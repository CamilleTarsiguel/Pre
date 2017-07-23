%%
allscen=[601 603 604 606 607 608 610 612];
allscen=[602,605,611,613,615,617,619];
allscen=600;
for s=allscen
    try
        scenario=s;
        [metrics2d, metrics3d, allens, stateInfo]=dcTracker(scenario);
%         startPT=stateInfo;
        save(sprintf('/home/amilan/research/projects/dctracking/out/s%04d.mat',s),'stateInfo','metrics*');
    catch err
         fprintf('%s\n',err.identifier);
    end
end

%%
for s=allscen
    try
        load(sprintf('/home/amilan/research/projects/dctracking/out/s%04d.mat',s));
        tracklets=convertToKITTI(stateInfo);
        writeLabels(tracklets,'./out',s-600);
    catch err
        fprintf('%s\n',err.identifier);
    end
end

%%
meanmets=[];
for s=allscen
    load(sprintf('/home/amilan/research/projects/dctracking/out/s%04d.mat',s));
    meanmets(end+1,:)=metrics2d;
    printMetrics(metrics2d,metricsInfo2d,0);
end
meanmets=mean(meanmets);
meanmets(4:11)=round(meanmets(4:11));
fprintf('\n');
printMetrics(meanmets);