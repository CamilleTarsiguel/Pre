function [metrics2d, metrics3d, allens, stateInfo, sceneInfo]= ...
    dcTracker(scene,options)
    
% Discrete-Continuous Optimization for Multi-Target Tracking
% with exclusion modeling
%
% This code served as basis for the following papers
%
% 
%  
%  Multi-Target Tracking by Discrete-Continuous Energy Minimization
%  A. Milan, K. Schindler and S. Roth. submitted to PAMI
%  
%  Detection- and Trajectory-Level Exclusion in Multiple Object Tracking
%  A. Milan, K. Schindler and S. Roth. In CVPR 2013 
%
%  Discrete-Continuous Optimization for Multi-Target Tracking
%  A. Andriyenko, K. Schindler and S. Roth. In CVPR 2012
%  
%
%  
% Please cite a subset of the publications if you use the code 
%  
%  
%  
% (C) Anton Milan, 2011-2015
%     anton.milan at adelaide dot edu dot au
%
% The code may be used free of charge for non-commercial and
% educational purposes, the only requirement is that this text is
% preserved within the derivative work. For any other purpose you
% must contact the authors for permission. This code may not be
% redistributed without written permission from the authors.
%  
%  
%  
% input:
% scene - description of the present scene, either
%               a numerical value representing a specific video sequence or
%               a struct containing all info or
%               a file name to the scene.ini
% 
% options - a struct with all necessary options and parameters
%
% output
% metrics2d - quantitative metrics computed as bounding box overlap 
% metrics3d - quantitative metrics computed as euclidean distance on ground plane
% alles     - energy values
% stateInfo - a struct with tracking results

% global start time
global dcStartTime
dcStartTime=tic;

% for compact output
format compact

% add all necessary paths
addPaths;

%% declare global variables
global detections nPoints sceneInfo opt globiter gtInfo scenario;
globiter=0;

global Dcost Lcost Scost

global TNeighbors SNeighbors T% TODO REMOVE
global alldpoints; % TODO REMOVE THIS LINE
global stInfTMP enTMP

global LOG_allens LOG_allmets2d LOG_allmets3d %for debug output
global LOG_allensdetailed
global tokeep mhs_% debugging
global mhs used labeling outlierLabel

% what dataset/sequence?
scenario=41;
if nargin, scenario=scene; end

% fill options struct with default if not given as parameter
% opt=getDCOptions;
if nargin>1    
    if isstruct(options)
        opt=options;
    elseif ischar(options)
        opt=readDCOptions(options);
%         opt
    else
        error('options parameters not recognized')
    end
else
    opt=readDCOptions('config/default2d.ini');    
end

randrun=opt.randrun;


%% seed random for deterministic results
if verLessThan('matlab','7.11'),    rand('seed',randrun);    randn('seed',randrun);
else    rng(randrun);
end


% fill scene info
sceneInfo=getSceneInfo(scenario);

% fill in empty output
[metrics2d, metrics3d]=getMetricsForEmptySolution();
stateInfo=getStateForEmptySolution(sceneInfo,opt);
allens=0;
LOG_allens=[];
LOG_allensdetailed=[];

% sceneInfo=getSceneInfoDCDemo;

frames=1:length(sceneInfo.frameNums);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frames=21:30; % do a part of the whole sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('doframes2.txt','file'), frl=load('doframes2.txt'); frames=frl(1):frl(2); end

if isfield(opt,'frames'), frames=opt.frames; end
if length(sceneInfo.frameNums)<frames(end), frames=1:length(sceneInfo.frameNums); end % if section outside range

%sceneInfo.frameNums=sceneInfo.frameNums(frames);




%% remove unnecessary frames from GT
%gtInfo=cropFramesFromGT(sceneInfo,gtInfo,frames,opt);

%
if opt.visOptim,  reopenFig('optimization'); end
%% load detections
[detections, nPoints]=parseDetections(sceneInfo,frames,opt.detThreshold);
[detections, nPoints]=cutDetections(detections,nPoints,sceneInfo, opt);
% detMatrices=getDetectionMatrices(detections);

% additional scene Info from detections
if ~isfield(sceneInfo,'targetAR')
    sceneInfo.targetAR=mean([detections(:).wd]./[detections(:).ht]);
end
if ~isfield(sceneInfo,'targetSize')
    sceneInfo.targetSize=mean([detections(:).wd]./2);
    if opt.track3d
        sceneInfo.targetSize=350; % 35 cm
    end
end

T=size(detections,2);                   % length of sequence
stateInfo.F=T; stateInfo.frameNums=sceneInfo.frameNums;

% fill remaining parameters that depend on sceneInfo
conffile='config/conOpt.txt';
global glclnr
if isdeployed
    conffile=sprintf('config/%d/conOpt.txt',glclnr);
end
opt=getAuxOpt(conffile,opt,sceneInfo,T);
checkDCOptions(opt);    % Check options for correctness

% degenerate case, sequence has less than 5 detection
% return empty solution
if numel([detections(:).xi])<5
    fprintf('Too few detections present. Exit.\n');
    [metrics2d, metrics3d, m2i, m3i, addInfo2d, addInfo3d]=getMetricsForEmptySolution();
%     stateInfo.Xi=[];stateInfo.Yi=[];stateInfo.Xgp=[];stateInfo.Ygp=[];stateInfo.X=[];stateInfo.Y=[];
%     stateInfo.X=[];stateInfo.Y=[];
    stateInfo.sceneInfo=sceneInfo;    stateInfo.opt=opt;
    stateInfo.splines=[];    stateInfo.outlierLabel=0;    stateInfo.labeling=[];
    
    allens=zeros(1,5);
    return;
end



% top image limit
sceneInfo.imTopLimit=max(1,min([detections(:).yi]-10));

% image plane frustum on the image plane
sceneInfo=computeImBordersOnGroundPlane(opt,sceneInfo,detections);


%% put all detections into a single vector
alldpoints=createAllDetPoints(detections);
nPoints=length(alldpoints.xp);  % number of detections |D|

%% print header
printHeader(sceneInfo,scenario,randrun);
fprintf('--- Discrete-Continuous Multi-target Tracking ---\n');
printSceneInfo;
printParams(opt);


%% create spatio-temporal neighborhood graph
TNeighbors=getTemporalNeighbors(alldpoints);
SNeighbors=getSpatialNeighbors(alldpoints);
Nhood.TNeighbors=TNeighbors; Nhood.SNeighbors=SNeighbors;
global glNhood; glNhood=Nhood; % ONLY for vis

%% init solution
% generate initial spline trajectories
mhs=getSplineProposals(alldpoints,opt.nInitModels,T);

%% get splines from EKF
for ekfexp=opt.startFromEKF
    dbfolder=fullfile('data','init'); 
    initsolfile=fullfile(dbfolder,'ekftracking',sceneInfo.dataset,sceneInfo.sequence,sprintf('e%04d.mat',ekfexp));
    %     initsolfile
    mhsekf=getSplinesFromEKF(initsolfile,frames,alldpoints,T);
    mhs=[mhs mhsekf];
end

global startPT
%% get splines from DP [Pirsiavash et al.]
if opt.startFromPir     
        pOpt=getPirOptions;
        startPT=runDP(detections,pOpt,opt);
    
        if opt.track3d
            [startPT.X,startPT.Y]=projectToGroundPlane(startPT.Xi,startPT.Yi,sceneInfo);
            startPT.Xgp=startPT.X;startPT.Ygp=startPT.Y;
        else
            if sceneInfo.yshift
                startPT.Y=startPT.Y-startPT.H/2;
            end
        end

        if opt.track3d && opt.cutToTA,        startPT=cutStateToTrackingArea(startPT);    end
%         startPT=cropFramesFromGT(sceneInfo,startPT,frames,opt);
                        
        
        mhsp=getSplinesFromGT(startPT.X,startPT.Y,frames,alldpoints,T);
        mhs=[mhs mhsp];
        
        %fprintf('Initial Result: \n');
        %printFinalEvaluation(startPT, gtInfo, sceneInfo, opt);
        
        clear startPT
end

    
%% remove unn. frames from camPar (UNUSED)
% if opt.track3d
%     if length(sceneInfo.camPar)>1, sceneInfo.camPar=sceneInfo.camPar(frames); end
% end

%% set initial labeling to all outliers
nCurModels=length(mhs);
nLabels=nCurModels+1; outlierLabel=nLabels;
labeling=nLabels*ones(1,nPoints); % all labeled as outliers
if opt.startFromGT
    labeling=getGreedyLabeling(alldpoints, gtSplines,T);
end


% unary is constant to outlierCost
Dcost=getUnarySpline(nLabels,nPoints,mhs,alldpoints,opt);   
Scost=opt.pairwiseFactor-opt.pairwiseFactor*eye(nLabels);
Lcost=getLabelCost(mhs,opt);
used=[]; % T*  is empy


% initial energy
energy=evaluateEnergy(alldpoints, Nhood, labeling, mhs, opt, []);

initialEnergy = energy;
bestEnergy = initialEnergy;
bestE=bestEnergy.value; E=energy.value;
initialEnergy=energyLogs(initialEnergy,0,0,0,toc(dcStartTime));
LOG_allens=[LOG_allens initialEnergy];
LOG_allmets2d=zeros(1,14);LOG_allmets3d=zeros(1,14);

% initial metrics
% if opt.startFromGT
%     [inMets2d, inMets3d]=printDCUpdateOGM(stateInfo,mhs,1:length(mhs),0,0,0,energy,'i');
% else
%     [inMets2d, inMets3d]=printDCUpdateOGM(stateInfo,mhs,[],0,0,0,energy,'i');
% end

%% first plot
drawDCUpdate(mhs,1:length(mhs),alldpoints,0,outlierLabel,TNeighbors,frames);

nAdded=0; nRemoved=0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% start energy minimization loop %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itcnt=0; % only count one discrete-continuous cycle as one iteration
iteachcnt=0; % count each discrete and each continuous optimization step
mhsafterrefit=mhs; labbeforeHSM=labeling;

m2d=zeros(1,14);m3d=m2d;
while 1
    oldN=length(mhs);
    
    % update info about when each spline was last used (active)
    for m=1:length(mhs)
        if ~isempty(intersect(m,used))
            mhs(m).lastused=0;
        else
            mhs(m).lastused=mhs(m).lastused+1;
        end
    end
    
    mhs_=mhs;
    tokeep=find([mhs.lastused]<opt.keepHistory);
    mhs=mhs(tokeep);
    
    nRemoved=oldN-length(tokeep);
    nCurModels=length(mhs); nLabels=nCurModels+1; outlierLabel=nLabels;
    
    
        
    % rearrange labeling, hacky --> improve TODO
    ul=setdiff(unique(labeling),length(mhs_)+1);
    for cl=ul
        newl=find(tokeep==cl);
        labeling(labeling==cl)=newl;
        used(used==cl)=newl;
    end
    labeling(labeling==length(mhs_)+1)=outlierLabel;
    
    energy_=energy;
    
    %% Discrete optimization - relabel all detections
    % minimize discrete Energy E(f), (Eq. ??)
    
    
    Dcost=getUnarySpline(nLabels,nPoints,mhs,alldpoints,opt);
    Lcost=getLabelCost(mhs,opt);
    Scost=opt.pairwiseFactor-opt.pairwiseFactor*eye(nLabels);
    precomp.Dcost=Dcost; precomp.Lcost=Lcost;
        
        
    % save energies before discrete minimization, type=1
    saveEnergySnapshot(alldpoints, detections, Nhood, labeling, mhs, opt, sceneInfo, 1);
    
        
    ticAlphaExp=tic;
    [~, labogm]=alphaExpansion(Dcost,Scost,Lcost, labeling,...
        Nhood, alldpoints, mhs,   opt);
    tocAlphaExp=toc(ticAlphaExp);
    
    for id=1:length(mhs)
        [mhs(id).labelCost, mhs(id).lcComponents] = ...
            getSplineGoodness(mhs(id),1,alldpoints,opt.seqLength);
    end
    Dcost=getUnarySpline(nLabels,nPoints,mhs,alldpoints,opt);
    Lcost=getLabelCost(mhs,opt);
    Scost=opt.pairwiseFactor-opt.pairwiseFactor*eye(nLabels);
    precomp.Dcost=Dcost; precomp.Lcost=Lcost;
    
    
    % retrieve new labeling
    labeling=    double(labogm);
    used=setdiff(unique(labeling),outlierLabel);
    
    % What is the energy after discrete optimization?    
    energy=evaluateEnergy(alldpoints, Nhood, labeling, mhs, opt, []);
    
    % count one iteration
    iteachcnt=iteachcnt+1;
    energy = energyLogs(energy, 1,iteachcnt,tocAlphaExp,toc(dcStartTime));
    LOG_allens=[LOG_allens energy];
    
    % print status
%     [m2d, m3d]=printDCUpdateOGM(stateInfo,mhs,used,nAdded,nRemoved,iteachcnt,energy,'d');
    %      pause
    
    
    %% Try removing trajectories (labels), one by one
    ticGreedy=tic;
    [labeling, used]=...
        greedyRemoval(alldpoints, mhs, labeling, outlierLabel, Nhood, opt, energy);
    tocGreedy=toc(ticGreedy);
    
    energy=evaluateEnergy(alldpoints, Nhood, labeling, mhs, opt, []);    
    E=energy.value;
    
    
    % if new energy worse, same or only slightly better,
    % or max iter reached, restore previous labeling and done
    if E+E*1e-3 >= bestE || iteachcnt >= opt.maxItCount
        if iteachcnt >= opt.maxItCount
            printMessage(2, 'Number of maximal iterations reached.\n');
        else
            printMessage(2, 'Discrete Optimization did not find a lower energy.\n');
        end
        
        mhs=mhsafterrefit;
        labeling=labbeforeHSM;
        
        energy=energy_;
        nCurModels=length(mhs); nLabels=nCurModels+1; outlierLabel=nLabels;
        used=setdiff(unique(labeling),outlierLabel);
        break;
    end
   
    
    % otherwise refit and adjust models
    bestE=energy.value;
    bestEnergy=energy;
    
    % count one iteration
    itcnt=itcnt+1;
    iteachcnt=iteachcnt+1;
    energy=energyLogs(energy,2,iteachcnt,tocGreedy,toc(dcStartTime));
    LOG_allens=[LOG_allens energy];
    
    % save energies before continuous minimization, type=2
    saveEnergySnapshot(alldpoints, detections, Nhood, labeling, mhs, opt, sceneInfo, 2);
    
    
    outlierLabel=nLabels;
    used=setdiff(unique(labeling),outlierLabel);
    
    
    % print update
    if opt.visOptim,  reopenFig('optimization'); end
    drawDCUpdate(mhs,used,alldpoints,labeling,outlierLabel,TNeighbors,frames);
    
    if opt.visOptim==2,  reopenFig('labeling');
        prepFigure; drawPoints(alldpoints,labeling,outlierLabel,TNeighbors);
        drawnow;
    end
    
    
    % print update
    [m2d, m3d]=printDCUpdateOGM(stateInfo,mhs,used,nAdded,nRemoved,iteachcnt,energy,'r');

        
    LOG_allmets2d(iteachcnt,:)=m2d;LOG_allmets3d(iteachcnt,:)=m3d;
    
    %% Continuous minimization
    % Rrefit trajectories (Eq. ??)
    mhsbeforerefit=mhs;
    
    ticCont=tic;
    stInfTMP=stateInfo; enTMP=energy; % TEMPORARY, Debugging
    mhsnew=minContinuousEnergy( ...
        mhs,opt, alldpoints, Nhood, used, labeling, Dcost, sceneInfo, stateInfo, energy);
    tocCont=toc(ticCont);
    
    
    %     mhsnew=reestimateSplines(allpoints,used,labeling,minCPs,ncpsPerFrame,1);
    mhsafterrefit=mhsnew;
    labbeforeHSM=labeling; % labeling before modifying hypotheses space
    
    Dcost=getUnarySpline(nLabels,nPoints,mhsnew,alldpoints,opt);
    Lcost = getLabelCost(mhsnew,opt);
    Scost=opt.pairwiseFactor-opt.pairwiseFactor*eye(nLabels);
    precomp.Dcost=Dcost; precomp.Lcost=Lcost;
    mhs(used)=mhsnew(used);
    nCurModels=length(mhs);
        
    energy=evaluateEnergy(alldpoints, Nhood, labeling, mhs, opt, []);
    
    
    % count one iteration
    iteachcnt=iteachcnt+1;
    
    energy=energyLogs(energy,3,iteachcnt,tocCont,toc(dcStartTime));
    LOG_allens=[LOG_allens energy];
    
    if energy.value>bestEnergy.value
        warning('Something wrong with continuous energy minimization. ...\nBest: %.2f, new: %.2f\n',bestE, energy.value);
    end
    bestE=energy.value;
    bestEnergy=energy;
    
    %     clear Scost Dcost Lcost
    
    % print update
    %     drawDCUpdate(mhs,1:length(mhs),alldpoints,0,outlierLabel,TNeighbors,frames);
    if opt.visOptim,  reopenFig('optimization'); end
    drawDCUpdate(mhs,used,alldpoints,labeling,outlierLabel,TNeighbors,frames);
    %     pause
    
    [m2d, m3d]=printDCUpdateOGM(stateInfo,mhs,used,nAdded,nRemoved,iteachcnt,energy,'c',1);
    %     LOG_allens(iteachcnt,:)=double([Eun Epw Eex Ela Eprox]);
    LOG_allmets2d(iteachcnt,:)=m2d;LOG_allmets3d(iteachcnt,:)=m3d;
    
    
    
    %% Expand the hypothesis space
    if nCurModels<opt.maxModels
        nModelsBeforeAdded=nCurModels;
        
        % merge crisscross
        mhsnew=mergeCrissCross(mhs,alldpoints,T,opt);
        mhs=[mhs mhsnew];
        
        %% get random new proposals
        mhsnew=getSplineProposals(alldpoints,opt.nAddRandomModels,T);
        mhs=[mhs mhsnew];
        
        %% get new proposals from outliers
        outlierPoints=find(labeling==outlierLabel); % indexes
        if length(outlierPoints)>4
            outlpts=selectPointsSubset(alldpoints,outlierPoints);
            mhsnew=getSplineProposals(outlpts,opt.nAddModelsFromOutliers,T);
            mhs=[mhs mhsnew];
        end
        
        %% extend existing
        mhs=extendSplines(alldpoints,mhs,used,labeling,T,E);
        
        %% shrink existing
        mhs=shrinkSplines(alldpoints,mhs,used,labeling,T,E);
        
        %% merge existing
        mhs=mergeSplines(alldpoints,mhs,used,labeling,T,E);
        
        %% split existing
        mhs=splitSplines(alldpoints,mhs,used,labeling,T,E);
        
        %% extend all to enable safety margin
        mhs=extendSplines2(alldpoints, mhs, used,opt);
    end
    nCurModels=length(mhs); nLabels=nCurModels+1; outlierLabel=nLabels;
    nAdded=nCurModels-length(mhsbeforerefit);
    
    if opt.visOptim==2,  reopenFig('hypotheses');
        prepFigure; drawSplines(mhs,1:length(mhs),0,[],frames)
        drawPoints(alldpoints,0,0,TNeighbors);
    end
    newl=labeling; newl(newl==nModelsBeforeAdded+1)=outlierLabel;
    labeling=newl;
    
end
% basically we are done
printMessage(1,'All done (%.2f min = %.2fh = %.2f sec per frame)\n',toc(dcStartTime)/60,toc(dcStartTime)/3600,toc(dcStartTime)/stateInfo.F);

iteachcnt=iteachcnt+1;

%% enumerate labeling to 1:N
newLabeling=zeros(1,nPoints);
for l=used, newLabeling(labeling==l)=find(used==l); end

% special handling of outlierLabel
lo=find(labeling==outlierLabel);
outlierLabel=length(used)+1;
if ~isempty(lo),    newLabeling(lo)=outlierLabel; end


%% final plot
% if usejava('desktop') && T<50
%     reopenFig('optimization');
%     prepFigure; drawPoints(alldpoints,newLabeling,outlierLabel,TNeighbors);
%     drawSplines(mhs(used),1:length(used),newLabeling,alldpoints,frames);
% end

%%
stateInfo=getStateFromSplines(mhs(used), stateInfo);
stateInfo=postProcessState(stateInfo);
if opt.remOcc,    stateInfo=removeOccluded(stateInfo); end
stateInfo.sceneInfo=sceneInfo;
stateInfo.opt=opt;
stateInfo.splines=mhs(used);
stateInfo.outlierLabel=outlierLabel;
stateInfo.labeling=newLabeling;
%% if we have ground truth, evaluate results
[metrics2d, metrics3d]=printFinalEvaluation(stateInfo, gtInfo, sceneInfo, opt);
% allens=double(LOG_allens(end,:));
% allens=0;
allens=[energy.value, energy.smoothness, ...
    energy.detExclusion, energy.labelCost.value, energy.trajExclusion];
    
%  fprintf('DO NOT FORGET LOGS\n');




% % %
% you can display the results with
% displayTrackingResult(sceneInfo,stateInfo)
