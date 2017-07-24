function [metrics, metricsInfo, additionalInfo, metrics3d, metricsInfo3d, additionalInfo3d]= ...
    evaluateDetections(detMatrices,gtInfo, sceneInfo, opt)
%

% global sceneInfo opt
if sceneInfo.gtAvailable
    detInfo=detMatrices; detInfo.frameNums=sceneInfo.frameNums;
    detInfo.X=detInfo.Xi; detInfo.Y=detInfo.Yi;
    [detInfo.F, detInfo.N]=size(detInfo);


%     detInfo
    printMessage(1,'\nDetections Evaluation (2D):\n');
%     gtInfo
%     detInfo
    [metrics, metricsInfo, additionalInfo]=CLEAR_MOT(gtInfo,detInfo,struct('eval3d',0));
    printMetrics(metrics,metricsInfo,1,[1 2 3 8 9]);
    printMessage(1,'\n');
    
    if opt.track3d
        if opt.cutToTA
            error('cutToTA not implemented');
        end
        [detMatrices.Xgp, detMatrices.Ygp]= ...
            projectToGroundPlane(detMatrices.Xi, detMatrices.Yi, sceneInfo);
        detInfo.Xgp=detMatrices.Xgp; detInfo.Ygp=detMatrices.Ygp;
    printMessage(1,'\nDetections Evaluation (3D):\n');
%     gtInfo
%     detInfo
    [metrics3d, metricsInfo3d, additionalInfo3d]=CLEAR_MOT(gtInfo,detInfo,struct('eval3d',3));
    printMetrics(metrics3d,metricsInfo3d,1,[1 2 3 8 9]);
    printMessage(1,'\n');
        
        %% TAKE CARE OF CROPPING PROPERLY!!!
%         detInfo=cutStateToTrackingArea(detInfo);
%         gtInfo=cutGTToTrackingArea(gtInfo, sceneInfo);
%         [detInfo.X, detInfo.Y]=projectToGroundPlane(detInfo.X,detInfo.Y,sceneInfo);
%         detInfo.Xgp=detInfo.X; detInfo.Ygp=detInfo.Y;
    end
    
end

end