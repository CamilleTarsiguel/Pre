function [labeling, used]= ...
    greedyRemoval(alldpoints, mhs, labeling, outlierLabel, Nhood, opt, energy)
% Remove a complete label one by one if it reduces the overall energy (Sec. 4.1, PAMI)

used=setdiff(unique(labeling),outlierLabel);

nLabels=length(mhs)+1; nPoints=length(alldpoints.xp);
Dcost=getUnarySpline(nLabels,nPoints,mhs,alldpoints,opt);
Lcost=getLabelCost(mhs,opt);
precomp.Dcost=Dcost;
precomp.Lcost=Lcost;

mcnt=0;
for labtorem=used(randperm(length(used)))
    mcnt=mcnt+1;
    trylab=labeling;
    trylab(trylab==labtorem)=outlierLabel;
    
    %         energy_r=evaluateFG(trylab,Dcost,Lcost, TNeighbors,opt.pairwiseFactor,SNeighbors,opt.exclusionFactor,tmpState);
    energy_r=evaluateEnergy(alldpoints, Nhood, trylab, mhs, opt, precomp);
    %         assert(isequal(energy_r,energy_r2));
    Eogmr=energy_r.value;
    %         tmpprox=proxcostGlob;
    %         tmpprox(mcnt,:)=0; tmpprox(:,mcnt)=0;
    %         Eproxr=sum(tmpprox(:));
    %         Eogmr=Eogmr+Eproxr;
    
    
    % if new energy is indeed lower ...
    if Eogmr<energy.value
        %             fprintf('before removal:%21.1f|%7.1f|%5.1f|%7.1f|%7.1f|%7.1f|%6.1f|%6.1f|%6.1f|%6.1f|%6.1f|\n', ...
        %                 energy.value, energy.data, energy.smoothness, energy.detExclusion, energy.trajExclusion, ...
        %                 energy.labelCost.value, energy.labelCost.regularizer, energy.labelCost.linVelocity,  ...
        %                 energy.labelCost.angVelocity, energy.labelCost.persistence, energy.labelCost.occlusionGaps);
        %
        fprintf('  %4i removed:%21.1f|%7.1f|%5.1f|%7.1f|%7.1f|%7.1f|%6.1f|%6.1f|%6.1f|%6.1f|%6.1f|\n', labtorem, ...
            energy_r.value, energy_r.data, energy_r.smoothness, energy_r.detExclusion, energy_r.trajExclusion, ...
            energy_r.labelCost.value, energy_r.labelCost.regularizer, energy_r.labelCost.linVelocity,  ...
            energy_r.labelCost.angVelocity, energy_r.labelCost.persistence, energy_r.labelCost.occlusionGaps);
        
        % ... set new energy value, labeling, and active trajectories
        energy=energy_r;
        labeling=trylab;
        used=setdiff(used,trylab);
    end
end

end