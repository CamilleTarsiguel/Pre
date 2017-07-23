function splines=refitSplines(splines,allpoints,Nhood, labeling, outlierLabel, opt)


global glclnr glexpnr % DEBUGGING

used=setdiff(unique(labeling),outlierLabel);
% reassign labeling to only use the active ones
nPoints=length(allpoints.xp);
newLabeling=zeros(1,nPoints);
for l=used, newLabeling(labeling==l)=find(used==l); end
newLabeling=labeling;

% [proxGlob, proxtGlob, proxcostGlob]=getSplineProximity(splines(lu),T,stInfo,opt);

allpointswithvpoints=generateVirtualPoints(allpoints, newLabeling,used);

splinesActive=splines(used);


for id=used
    [splines(id).labelCost, splines(id).lcComponents] = ...
        getSplineGoodness(splines(id),1,allpoints,opt.seqLength);
end
precomp.Lcost = getLabelCost(splines,opt);

tmpState.F=opt.seqLength;
tmpState=getStateFromSplines(splines(used), tmpState, ~opt.track3d);
[~, ~, proxcostGlob]=getSplineProximity(splines(used),opt.seqLength, tmpState, opt);
precomp.proxcost = proxcostGlob;


energy=evaluateEnergy(allpoints, Nhood, labeling, splines, opt, precomp);
mcnt=0;
for id=used
    mcnt=mcnt+1;
    
    oldSpline=splines(id);
    oldLcost=precomp.Lcost;
    oldPcost=precomp.proxcost;
    
    supportPts=find(newLabeling==id);
    xy=[allpointswithvpoints.xp(supportPts);allpointswithvpoints.yp(supportPts)];
    tr=allpointswithvpoints.tp(supportPts);
    confs=allpointswithvpoints.sp(supportPts);
	if length(unique(tr)) < 4
		continue;
	end
	
    
%     energy=evaluateEnergy(allpoints, Nhood, labeling, splines, opt, precomp);
    
%     save(sprintf('config/%d/sfit-debug-%d',glclnr,glexpnr),'*');
    tryfit=splinefit(tr,xy,oldSpline.pieces,4,confs,'s');
    
    % singular matrix case
    if any(isnan(tryfit.coefs(:)))
        continue;
    end
    tryfit=adjustSplineStruct(tryfit, min(tr), max(tr), allpoints, opt.seqLength, 0, [], [], []);
%         tryfit=adjustSplineStruct(tryfit, oldSpline.start, oldSpline.end, allpoints, opt.seqLength, 0, [], [], []);
    
    
    splines(id)=tryfit;
    
    [prox2, proxt2, proxcost2]= ...
        getOneSplineProximity(splines(used),opt.seqLength, [], mcnt,opt);
    
    
    proxcostNew=proxcostGlob;
    if mcnt==1
        proxcostNew(1,2:end)=proxcost2(2:end);
    elseif mcnt==length(used)
        proxcostNew(1:end-1,mcnt)=proxcost2(1:end-1)';
    else
        proxcostNew(1:mcnt-1,mcnt)=proxcost2(1:mcnt-1);
        proxcostNew(mcnt,mcnt+1:end)=proxcost2(mcnt+1:end);
    end
    precomp.proxcost = proxcostNew;
    
    [splines(id).labelCost, splines(id).lcComponents] = ...
        getSplineGoodness(splines(id),1,allpoints,opt.seqLength);
    precomp.Lcost = getLabelCost(splines,opt);
    
    
%     energyNew_=evaluateEnergy(allpoints, Nhood, labeling, splines, opt, []);
    energyNew=evaluateEnergy(allpoints, Nhood, labeling, splines, opt, precomp);
%     assert(isequal(energyNew,energyNew_));
    
    
    if energy.value < energyNew.value
        splines(id)=oldSpline;
        precomp.Lcost=oldLcost;
        precomp.proxcost=oldPcost;
    else
        proxcostGlob=proxcostNew;
        energy=energyNew;
    end

    
end

% splinesNewActive=splinesActive;
for id=used
    [splines(id).labelCost, splines(id).lcComponents] = ...
        getSplineGoodness(splines(id),1,allpoints,opt.seqLength);
end

% splines(used)=splinesNewActive;
