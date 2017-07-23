function energy = ...
    evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, precomp)

% Evaluate entire energy
%
% Graphical Model entities
% alldpoints - detection points
% Nhood.TNeighbors - Temporal (submodular) connections
% Nhood.SNeighbors - Spatial (non-submodular) connections
% labeling - current labeling

% Continuous entities
% splines - trajectories

% Shared entitites
% opt - weights
% T - sequence length

% precomp - a struct containing precomputed quantites



% first precompute all needed auxiliary quantities
nCurModels=length(splines); nLabels=nCurModels+1; outlierLabel=nLabels;
nPoints = length(alldpoints.xp);

b=labelcount(labeling,nLabels);
labelUsed=~~b;
labelUsed=double(labelUsed);
lu= find(labelUsed(1:end-1));

% if isempty(precomp)
    if ~isfield(precomp,'Dcost')
%         pause
        precomp.Dcost = ...
            getUnarySpline(nLabels,nPoints,splines,alldpoints,opt);
    end
    
    if ~isfield(precomp,'Lcost')
        for id=1:length(splines)
            [splines(id).labelCost, splines(id).lcComponents] = ...
                getSplineGoodness(splines(id),1,alldpoints,opt.seqLength);
        end
        precomp.Lcost = getLabelCost(splines,opt);
%         sum([splines(lu).labelCost])
%         sum(precomp.Lcost(lu))
        
%         for id=lu
%             [splines(id).labelCost, splines(id).lcComponents] = ...
%                 getSplineGoodness(splines(id),1,alldpoints,opt.seqLength);
%         end
%         precomp.Lcost = getLabelCost(splines(lu),opt);
    end    
    
    if ~isfield(precomp,'proxcost')
        tmpState.F=opt.seqLength;
        tmpState=getStateFromSplines(splines(lu), tmpState, ~opt.track3d);
%         lu
%         save('tmpvars.mat','splines','lu','opt','tmpState');
        [~, ~, proxcost]=getSplineProximity(splines(lu),opt.seqLength, tmpState, opt);
        precomp.proxcost = proxcost;
    end
% end
% precomp
% isfield(precomp,'proxcost')

UCost=precomp.Dcost;
LCost=precomp.Lcost;

PWCost=precomp.proxcost;
if size(precomp.proxcost,1) > length(lu)
    PWCost=precomp.proxcost(lu',:);
    PWCost=PWCost(:,lu);
end
% PWCost

% check precomputed quantites
assert(nLabels==size(UCost,1));
assert(nPoints==size(UCost,2));

%% unaries    
tmp=cumsum(nLabels*ones(1,nPoints))-nLabels;
labinds=tmp+labeling;
Eun=sum(UCost(labinds));

%% Potts for temporal smoothness
[pwu, pwv]=find(Nhood.TNeighbors);
Epw=opt.pairwiseFactor*sum(labeling(pwu)~=labeling(pwv));

%% Detection-Level Exclusion
[exu, exv]=find(Nhood.SNeighbors);
Eex=0;
if opt.exclusionFactor
    Eex=opt.exclusionFactor*sum(labeling(exu)==labeling(exv) & labeling(exu)~=outlierLabel);
end

%% Labelcost
Ela=xdotyMex(labelUsed,LCost);
% Ela=sum(LCost);


%% Trajectory-level exclusion
Eprox=sum(PWCost(:));

% % % global gpwcost
% % % gpwcost=PWCost;
% % % % splines=mhs(lu);
% % % activeSplines=splines(lu);
% % % EexcValue=0;
% % % if ~isempty(lu)
% % %     stateVec=splinesToVec(activeSplines);
% % %     if opt.mex
% % % 
% % %         parEexc=opt.conOpt.enParEexc; % Eexc par
% % %         spInfo=[
% % %             [activeSplines(:).pieces];
% % %             [activeSplines(:).start];
% % %             [activeSplines(:).end]]';   
% % %         br=[activeSplines(:).breaks];
% % % 
% % %         brInfo=ones(length(activeSplines),2);
% % %         for id=1:length(activeSplines)
% % %             brInfo(id,1)=activeSplines(id).breaks(1);
% % %             brInfo(id,2)=activeSplines(id).breaks(end);
% % %         end
% % %         Nsp=length(activeSplines);
% % %         pcindex=[];
% % %         brindex=[]; brindxcnt=1;
% % %         for id=1:Nsp
% % %             pcindex=[pcindex id*ones(1,activeSplines(id).pieces)];
% % %             brindex=[brindex brindxcnt:brindxcnt+activeSplines(id).pieces-1];
% % %             brindxcnt=brindxcnt+activeSplines(id).pieces+1;
% % %         end
% % %         indexexc=[pcindex; brindex];
% % %         indexexc=indexexc-1;
% % % 
% % %         [EexcValue, ~]=Eexc_mex(stateVec,parEexc,spInfo,brInfo,br,indexexc);
% % %     else
% % %         EexcValue=Eexc(stateVec,activeSplines);
% % % 
% % %     end
% % % end
% % % EexcValue = opt.proxcostFactor * EexcValue;
% % % Eprox2=EexcValue;
% % % 
% % % assert(abs(Eprox-Eprox2)<1,'hmm... %f %f %f\n',Eprox,Eprox2,Eprox-Eprox2)

%% combine all
activeSplines=splines(lu);
allcomponents=[activeSplines.lcComponents];
allcomponents=reshape(allcomponents,5,length(lu))';
components=sum(allcomponents,1);

E=Eun+Ela+Epw+Eex+Eprox;

energy.value = E;
energy.data = Eun;
energy.smoothness = Epw;
energy.detExclusion = Eex;
energy.trajExclusion = Eprox;
energy.labelCost.value = Ela;

energy.labelCost.regularizer = opt.labelCost*length(lu);
energy.labelCost.linVelocity = components(1);
energy.labelCost.angVelocity = components(2);
energy.labelCost.persistence = components(3);
energy.labelCost.occlusionGaps = components(4);
energy.labelCost.segConsistency = components(5);


