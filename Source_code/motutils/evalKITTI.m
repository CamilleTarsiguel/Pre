function metricsKITTI=evalKITTI(infos,setting,cls)
    % evaluate KITTI Training Set
    % requires a struct array infos 1x21
    
    metricsKITTI=zeros(1,22);
    if length(infos)~=21
        fprintf('wrong data');
        return;
    end
    
    if nargin<2, setting='a'; end
    if nargin<3, cls='car'; end
    
    cls=lower(cls);
    resdir=sprintf('results/%s/data',setting);
    
    thiswd=pwd;
    
    pathToKittiDevkit='../motutils/external/KITTI/devkit_tracking/python';
    cd(pathToKittiDevkit)
    if ~exist(resdir,'dir'), mkdir(resdir); end
    
    allscen=1:21;
    for scenario=allscen
        stateInfo=infos(scenario).stateInfo;
        cd(thiswd)
        tracklets=convertToKITTI(stateInfo,cls);
        cd(pathToKittiDevkit)
	addpath(genpath('..'));

        writeLabels(tracklets,resdir,scenario-1);
    end
    
    pythonstr=sprintf('!~/software/python/bin/python evaluate_tracking.py %s',setting)
    eval(pythonstr);
    metricsKITTI=dlmread(sprintf('results/%s/stats_%s.txt',setting,cls));
    
    cd(thiswd);
end
