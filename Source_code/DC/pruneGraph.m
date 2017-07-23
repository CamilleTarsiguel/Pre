function [pPoints,pTN, pSN, pDcost, pLabeling, keeppts]= ...
    pruneGraph(allpoints,labeling,Dcost,TNeighbors,SNeighbors,tmpState,expLab,outlierLabel,mh,opt,T)
% discard portions of graph that are irrelevant for current move


allpind=1:length(allpoints.xp);
keeppts=allpind;
pPoints=allpoints;
pTN=TNeighbors;
pSN=SNeighbors;
pDcost=Dcost;
pLabeling=labeling;

% do not prune if we are expanding on the outlier label
if expLab==outlierLabel
    return;
end



X=tmpState.X(:,expLab);Y=tmpState.Y(:,expLab);

% include points beyond temporal life span
exthead=2;
exttail=2;
splstart=mh.start-exttail;
splend=mh.end+exthead;

% which points are in temporal life span
inSplineTimespan=find(allpoints.tp>=splstart & allpoints.tp<=splend);

pt=allpoints.tp(inSplineTimespan);
xt=[X(pt) Y(pt)]; %pts on spline   
px=allpoints.xp(inSplineTimespan);
py=allpoints.yp(inSplineTimespan);


datapts=[px; py]';

% compute distances on spline to points
alldists=xt-datapts;
    
alldists=alldists';
alldists=sqrt(sum(alldists.^2)); % L2 norm in meter or pixels


% keep those close to spline
keepptsTS=find(alldists<opt.tau*5);
keepptsAbs=inSplineTimespan(keepptsTS);
keeppts=keepptsAbs;

pPoints.xp=allpoints.xp(keepptsAbs);pPoints.yp=allpoints.yp(keepptsAbs);
pPoints.sp=allpoints.sp(keepptsAbs);pPoints.tp=allpoints.tp(keepptsAbs);


% adjust neighborhoods
pTN=TNeighbors(keepptsAbs',:);
pTN=pTN(:,keepptsAbs);

pSN=SNeighbors(keepptsAbs',:);
pSN=pSN(:,keepptsAbs);

% adjust unary costs
pDcost=Dcost(:,keepptsAbs);

% ... and label costs
pLabeling=labeling(keepptsAbs);



end