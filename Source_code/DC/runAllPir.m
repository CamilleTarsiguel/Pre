%%
allscen=[500:520 550:578];
% allscen=603;
% allscen=550;
allscen=600:620;
allscen=600;
% allscen=508;
for s=allscen
    try
        scenario=s;
        [metrics2d, metrics3d, allene, stateInfo]=test_pirsiavash(s,pOpt,opt);
        startPT=stateInfo;
        save(sprintf('/home/amilan/research/projects/dctracking/data/init/dptracking/startPT-pir-s%04d.mat',s),'startPT');
    catch err
        fprintf('oh\n');
        fprintf('%s\n',err.identifier)
        pause
    end
end