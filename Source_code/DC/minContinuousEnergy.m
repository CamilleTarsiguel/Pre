function splinesnew=minContinuousEnergy( ...
    splines,opt, alldpoints, Nhood, used, labeling, Dcost, sceneInfo, stateInfo, energy)
% This function calls an iterative continuous optimization algorithm to minimizeIterativeBS
% the continuous part of the energy

[nLabels, nPoints]=size(Dcost);
outlierLabel=nLabels;

global LOG_allens dcStartTime

 energy1=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);
% splines=readjustSplineRange(splines,labeling,used, alldpoints,Nhood,opt);
%  energy2=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);
% energy1
% energy2
 
if strcmpi(opt.conOpt.alg,'lsfit')
%     splinesnew=reestimateSplines(alldpoints,Nhood,used,labeling,nLabels,splines,Dcost,opt,sceneInfo,1);
    splinesnew=minContinuousEnergySimple(splines, used, alldpoints, labeling, Nhood, outlierLabel, opt);
    energy3=evaluateEnergy(alldpoints, Nhood, labeling, splinesnew, opt, []);
% elseif strcmpi(opt.conOpt.alg,'simplex')
%     energy=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);
%     splinesnew1=reestimateSplines(alldpoints,Nhood,used,labeling,nLabels,splines,Dcost,opt,sceneInfo,0);
%     splinesnew=minContinuousEnergySimple(splines, used, alldpoints, labeling, Nhood, outlierLabel, opt);
%      energy1=evaluateEnergy(alldpoints, Nhood, labeling, splinesnew1, opt, []);
%       energy2=evaluateEnergy(alldpoints, Nhood, labeling, splinesnew2, opt, []);
%       splinesnew=splinesnew2;
%       [energy.value energy1.value energy2.value]

%       pause
else
%     energy=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);

%     if opt.readjustSE
%         ticRSR=tic;
%         splines=readjustSplineRange(splines,labeling,used, alldpoints,Nhood,opt);
%         tocRSR=toc(ticRSR);
%         energy1=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);
%         [m2d, m3d]=printDCUpdateOGM(stateInfo,splines,used,0,0,0,energy1,'a');
%         energy1=energyLogs(energy1,4,0,tocRSR,toc(dcStartTime));
%         LOG_allens=[LOG_allens energy1];
%     end
    
    energy2.value=0;
    if opt.conOpt.initSimple
        ticRefit=tic;
        splines=refitSplines(splines,alldpoints, Nhood,labeling, nLabels, opt);
        tocRefit=toc(ticRefit);        
        
        energy2=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);
        [m2d, m3d]=printDCUpdateOGM(stateInfo,splines,used,0,0,0,energy2,'s');
        energy2=energyLogs(energy2,5,0,tocRefit,toc(dcStartTime));
        LOG_allens=[LOG_allens energy2];

%             splines=minContinuousEnergySimple(splines, used, alldpoints, labeling, Nhood, outlierLabel, opt);
    end
    splinesnew = minimizeIterativeBS(alldpoints, splines, used, labeling, outlierLabel, nPoints, energy, opt, sceneInfo);
    energy3=evaluateEnergy(alldpoints, Nhood, labeling, splinesnew, opt, []);
%     [energy.value energy1.value energy2.value energy3.value]
end

% a quick check, if energy did not increase
% this is due to our recent switch from pp to b-spline representation
% within matlab
if energy3.value>=energy1.value
    splinesnew=splines;
end


end