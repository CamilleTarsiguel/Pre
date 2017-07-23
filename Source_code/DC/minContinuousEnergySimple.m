function splines=minContinuousEnergySimple( ...
    splines, used, alldpoints, labeling, Nhood, outlierLabel, opt)
% minimize simple energy without labelcost

if isempty(used)
    return;
end


% what is the outlier label? shift to -1
% outlierLabel=length(splines)+1;
% labeling(labeling==outlierLabel)=-1;

% splinesActive=splines(used);
% reassign labeling to only use the active ones
% newLabeling=zeros(1,nPoints);
% for l=used, newLabeling(labeling==l)=find(used==l); end

% pointsAndVPoints=generateVirtualPoints(alldpoints,newLabeling);
[pointsAndVPoints, newLabeling]=generateVirtualPoints(alldpoints,labeling, used);

fitnorm='a';
if strcmpi(opt.conOpt.alg,'lsfit'), fitnorm='s'; end

for id=used
    splineOld=splines(id);
    supportPts=find(newLabeling==id);
    
    xy=[pointsAndVPoints.xp(supportPts);pointsAndVPoints.yp(supportPts)];
    t=pointsAndVPoints.tp(supportPts);
    
    origSupportPts=find(labeling==id);
    torig=alldpoints.tp(origSupportPts);
    
    confs=pointsAndVPoints.sp(supportPts);
    
    uniqueT=unique(t);
    nuniqueT=length(uniqueT);
    


    energy=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);
%     id
%     energy.value
    if nuniqueT<4
        continue;
    else
        [t, sortedind]=sort(t);
        xy=xy(:,sortedind); confs=confs(sortedind);

        sfit=splinefit(t,xy,splineOld.pieces,splineOld.order,confs,fitnorm);
        sfit=adjustSplineStruct(sfit, min(torig), max(torig), alldpoints, opt.seqLength, 0, [], [], []);        

    end
    
%     if id==43
%         xy
%         t
%         confs
%         sfit
%         sfit.coefs
%     end
    splines(id)=sfit;
    energyNew=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);
%     energy.value
%    pause
    if energy.value<energyNew.value
        splines(id)=splineOld;
    else
%         energy
%         energyNew
%         fprintf('Number %d refit, old: %f, new %f, diff: %f\n',id,energy.value,energyNew.value, energy.value-energyNew.value)
    end
    
end