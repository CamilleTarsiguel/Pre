%% view KITTI test results
allscen=[500:520 550:578];
allscen=[550:578];
for s=allscen
    fprintf('Sequence: %d\n',s);
    try
        load(sprintf('out/s%04d.mat',s));
        displayTrackingResult(stateInfo.sceneInfo,stateInfo);
    catch err
        fprintf('%s\n',err.identifier);
    end
end