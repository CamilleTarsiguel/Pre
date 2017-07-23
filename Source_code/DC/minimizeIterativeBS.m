function splinesNew = ...
    minimizeIterativeBS(alldpoints, splines, used, labeling, outlierLabel, nPoints, energy, dcOpt, sceneInfo)
% Perform continuous energy minimization (B-splines)
% Note that Smoothness, detection-level exclusion
% and hreg remain constant.

newOrigLabeling=zeros(1,nPoints);
% reassign labeling to only use the active ones
for l=used, newOrigLabeling(labeling==l)=find(used==l); end
newOrigLabeling(~newOrigLabeling)=-1;

global ptsAndVPoints
[ptsAndVPoints, labelingV]=generateVirtualPoints(alldpoints, labeling, used);
% ptsAndVPoints=alldpoints; labelingV=labeling;

global splinesActive splinesNewActive
global oldEnergy
oldEnergy=energy;

global T
T=dcOpt.seqLength;






%% precompute some values to speed up
% for id=used    
%     idPts=labeling==id;
%     torig=alldpoints.tp(idPts);
%     splines(id).start=min(torig);
%     splines(id).end=max(torig);
% end

for id=1:length(splines)
    spl=splines(id).bspline;
    spl.start=splines(id).start;
    spl.end=splines(id).end;
    bspl(id)=spl;
%     s1=splinesActive(id).start;
%     e1=splinesActive(id).end;
%     xypp=ppval(splinesActive(id),s1:e1);
%     xysp=spval(bspl(id),s1:e1);
%     sum(sum(xypp-xysp))
end



for id=1:length(splines)
    spl=splines(id);
%     idPts=labeling==id;
%     torig=alldpoints.tp(idPts);
%     tr=min(torig):max(torig);
    tr=spl.start:spl.end;
    
    % save piece indexes into splines structs
%     [~, index] = histc(tr,[-inf,spl.breaks(2:spl.pieces),inf]);
    [~, index] = histc(tr,[-inf,spl.bspline.knots(spl.bspline.order+1:spl.bspline.number),inf]);
    splines(id).index=index;
    
    % also save detections indices
    % which detection belongs to which trajectory
    % ORIG POINTS
    idPts= labeling==id;    
    alldpts.tp=alldpoints.tp(idPts);
    tr=alldpts.tp;
%     [~, index] = histc(tr,[-inf,spl.breaks(2:spl.pieces),inf]);
    [~, index] = histc(tr,[-inf,spl.bspline.knots(spl.bspline.order+1:spl.bspline.number),inf]);
    
%     if ~isequal(index,indexbspl)
%         error('spline index not right');
%     end
    
    splines(id).ptindex=index;

    % INCLUDING VPOINTS
    idPts= labelingV==id;    
    alldpts.tp=ptsAndVPoints.tp(idPts);
    tr=alldpts.tp;
    [~, index] = histc(tr,[-inf,spl.breaks(2:spl.pieces),inf]);
    splines(id).ptvindex=index;
    
    bspl(id).index=splines(id).index;
    bspl(id).ptindex=splines(id).ptindex;
end


% we only need active labels
global newLabeling
newLabeling=zeros(1,length(ptsAndVPoints.xp));
% reassign labeling to only use the active ones
for l=used, newLabeling(labelingV==l)=find(used==l); end

% what is the outlier label? shift to -1
% outlierLabel=length(splines)+1;
% labeling(labeling==outlierLabel)=-1;
newLabeling(newLabeling==0)=-1;


splinesActive=splines(used);



% mex stuff
global allt allos pts ptsv ptInfo
if dcOpt.mex
    [allt, allos, allos2, allkos]=detToIndexBS(splinesActive,newLabeling, ptsAndVPoints);
    ptsv=[ptsAndVPoints.xp; ptsAndVPoints.yp; ptsAndVPoints.sp; ptsAndVPoints.tp; newLabeling; allt; allos; allos2; allkos];
    

    [allt, allos, allos2, allkos]=detToIndexBS(splinesActive,newOrigLabeling, alldpoints);
    pts=[alldpoints.xp; alldpoints.yp; alldpoints.sp; alldpoints.tp; newOrigLabeling; allt; allos; allos2; allkos];
    
    ptInfo=zeros(T,2);
    for t=1:T
        tPts=find(alldpoints.tp==t);
        if ~isempty(tPts)
            ptInfo(t,:)=[tPts(1) tPts(end)];
        else
            ptInfo(t,:)=[0 -1];
        end
    end
    ptInfo=ptInfo-1;
end
% length([splines(:).ptindex])
% size(ptsv)
% pause
% add virtual points to fit on safety margin
% global oldAlldpoints;
% oldAlldpoints=alldpoints;
% global ptsAndVPoints

% size(alldpoints.xp)
% size(ptsAndVPoints.xp)
% size(labeling)
% size(newLabeling)
% pause



% xvec=splinesToVec(splinesActive);
% convert to B-form
% for id=1:length(splinesActive)
%     bspl(id)=mypp2sp(splinesActive(id));
% %     s1=splinesActive(id).start;
% %     e1=splinesActive(id).end;
% %     xypp=ppval(splinesActive(id),s1:e1);
% %     xysp=spval(bspl(id),s1:e1);
% %     sum(sum(xypp-xysp))
% end
% pause
bspl=bspl(used);
xvec=BSplinesToVec(bspl);


global gxvec
gxvec=xvec;
% global bspl newOrigLabeling
% [fold, df]=EconBS(xvec,splinesActive);
% checkgrad('EconBS',xvec, 1,bspl, ...
%         alldpoints, ptsAndVPoints, dcOpt, sceneInfo, newOrigLabeling, newLabeling, energy)

if strcmpi(dcOpt.conOpt.alg,'CGD')
%     save('evnars.mat','*');
    [X, ~, ~]= ...
        minimizeSplineBS(xvec,'EconBS', dcOpt.conOpt.maxIter,bspl, ...
        alldpoints, ptsAndVPoints, dcOpt, sceneInfo, newOrigLabeling, newLabeling, energy);
elseif strcmpi(dcOpt.conOpt.alg,'simplex')
    anonf = @(xvec)EconBS(xvec,bspl,alldpoints, ptsAndVPoints, dcOpt, sceneInfo, newOrigLabeling, newLabeling, energy);    
    options=optimset('MaxIter',dcOpt.conOpt.maxIter,'MaxFunEvals',100,'Display','off');
    [X] = fminsearch(anonf,xvec,options);    
elseif strcmpi(dcOpt.conOpt.alg,'fminunc')
    anonf = @(xvec)EconBS(xvec,bspl,alldpoints, ptsAndVPoints, dcOpt, sceneInfo, newOrigLabeling, newLabeling, energy);
    options = optimoptions('fminunc','MaxIter',dcOpt.conOpt.maxIter,'MaxFunEvals',100,'GradObj','on','Display','off');
    [X] = fminunc(anonf,xvec,options);
elseif strcmpi(dcOpt.conOpt.alg,'LBFGS')
    anonf = @(xvec)EconBS(xvec,bspl,alldpoints, ptsAndVPoints, dcOpt, sceneInfo, newOrigLabeling, newLabeling, energy);
    LBFGSoptimoptions = ...
        struct('GradObj','on','Display','off','MaxIter',dcOpt.conOpt.maxIter,'MaxFunEvals',400, ...
        'LargeScale','off','HessUpdate','bfgs','InitialHessType','identity','GoalsExactAchieve',1,'GradConstr',false);    
    [X] = fminlbfgs(anonf,xvec,LBFGSoptimoptions);
else
    warning('Unkown continuous minimizer specified. Performing fminunc');
    anonf = @(xvec)EconBS(xvec,bspl,alldpoints, ptsAndVPoints, dcOpt, sceneInfo, newOrigLabeling, newLabeling, energy);
    options=optimoptions('fminunc',dcOpt.conOpt.optimoptions);
    [X] = fminunc(anonf,xvec,options);
end



newxvec=X;
% [f, df] = EconBS(newxvec,splinesActive);
% fprintf('Start: %.2f,  final: %.2f\n',fold,f);

% remove virtual points
% alldpoints=oldAlldpoints;

splinesNewActive=splinesActive;

bsplNew=vecToBSplines(newxvec,bspl);
for id=1:length(bsplNew)
    bspl=bsplNew(id);
    pp=splinefit(bspl.start:bspl.end,spval(bspl,bspl.start:bspl.end),splinesActive(id).breaks,splinesActive(id).order);
%     bsplNew(id)
%     pp
%     splinesActive(id)
    assert(isequal(size(pp.coefs),size(splinesActive(id).coefs)));
    assert(isequal(pp.breaks,splinesActive(id).breaks));
    splinesNewActive(id).coefs=pp.coefs;
end


for id=1:length(splinesNewActive)
    [splinesNewActive(id).labelCost, splinesNewActive(id).lcComponents] = ...
        getSplineGoodness(splinesNewActive(id),1,alldpoints,T);
end

splinesNew=splines;
splinesNew(used)=splinesNewActive;

end