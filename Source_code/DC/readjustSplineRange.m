function splines=readjustSplineRange(splines,labeling,used, alldpoints,Nhood,opt)
% refit start and end points to labeled detections


for id=used
    [splines(id).labelCost, splines(id).lcComponents] = ...
        getSplineGoodness(splines(id),1,alldpoints,opt.seqLength);
end
precomp.Lcost = getLabelCost(splines,opt);

tmpState.F=opt.seqLength;
tmpState=getStateFromSplines(splines(used), tmpState, ~opt.track3d);
[~, ~, proxcostGlob]=getSplineProximity(splines(used),opt.seqLength, tmpState, opt);
precomp.proxcost = proxcostGlob;
global energy1 energy1_

energy=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, precomp);

mcnt=0;        
for id=used
    mcnt=mcnt+1;
    oldSpline=splines(id);
    oldLcost=precomp.Lcost;
    oldPcost=precomp.proxcost;
    
    
    idPts=labeling==id;
    t=alldpoints.tp(idPts);
    mint=min(t); maxt=max(t);
    
    if oldSpline.start > mint
%         energy=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, precomp);
        splines(id).start=mint;
        
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
            getSplineGoodness(splines(id),1,alldpoints,opt.seqLength);
        precomp.Lcost = getLabelCost(splines,opt);
    
    
%         energy1_=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);
        energy1=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, precomp);
%         assert(isequal(energy1_,energy1));

        
        % if new energy worse, restore old state
        if energy.value < energy1.value
             splines(id)=oldSpline;
            %             energy1=energy;
            precomp.Lcost=oldLcost;
            precomp.proxcost=oldPcost;
        else
            proxcostGlob=proxcostNew;
            energy=energy1;
        end
%             [id oldSpline.start mint]
%             energy
%             energy1
%             energy.labelCost
%             energy1.labelCost        
    end
    
    oldSpline=splines(id);
    oldLcost=precomp.Lcost;
    oldPcost=precomp.proxcost;
    
    if oldSpline.end < maxt
%         energy=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, precomp);
        splines(id).end=maxt;


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
            getSplineGoodness(splines(id),1,alldpoints,opt.seqLength);
        
        precomp.Lcost = getLabelCost(splines,opt);
        
        energy1=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, precomp);
%         energy1_=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);
%         assert(isequal(energy1,energy1_));
        
        if energy.value < energy1.value
             splines(id)=oldSpline;
             precomp.Lcost=oldLcost;
             precomp.proxcost=oldPcost;
        else
            proxcostGlob=proxcostNew;
            energy=energy1;
        end
%             [id oldSpline.end maxt]
%             energy
%             energy1
%             energy.labelCost
%             energy1.labelCost        
        
    end
%     pause
    
end