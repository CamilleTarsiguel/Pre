function paramAnal(jobname,jobid)

%% determine paths for config, logs, etc...
addpath(genpath('../motutils/'));
format compact

[~,hname]=system('hostname')
settingsDir=strsplit(jobname,'-');
runname=char(settingsDir{1})
learniter=char(settingsDir{2})
jid=char(settingsDir{3}); % unused

settingsDir=[runname '-' learniter];

confdir=sprintf('config/%s',settingsDir);

jobid=str2double(jobid);
confdir

resdir=sprintf('results/%s',settingsDir);
if ~exist(resdir,'dir'), mkdir(resdir); end
resdir
trainStartTime=tic;

resultsfile=sprintf('%s/res_%03d.mat',resdir,jobid);


try
    maxexper=dlmread(fullfile(confdir,'maxexper.txt'));
catch
    maxexper=20;
end

% if computed alread, just load it
if exist(resultsfile,'file')
  load(resultsfile);
else

  conffile=fullfile(confdir,sprintf('%04d.ini',jobid));
  conffile

  inifile=fullfile(confdir,'default.ini');
  inifile
  if ~exist(inifile,'file')
	  error('You must provide initial options file default.ini');
  end

  opt=readDCOptions(inifile);

  %% zip up logs (previous)
  if jobid==1
    prevSetting=sprintf('%s-%d',runname,str2double(learniter)-1)
    zipstr=sprintf('!sh ziplogs.sh %s',prevSetting)
    eval(zipstr);
  end


  % take care of parameters
  jobid

  rng(jobid);

  ini=IniConfig();
  ini.ReadFile(inifile);

  params=[];
  % we are only interested in [Parameters]
  sec='Parameters';
  keys = ini.GetKeys(sec);

  for k=1:length(keys)
      key=char(keys{k});
      params = [params ini.GetValues(sec,key)];
  end

  % param search dimension
  ndim=length(keys);
  lval=0; uval=2;

  if ndim==1
      pmult=uval*(jobid-1)/(maxexper-1);
      params=pmult*params;
  elseif ndim==2
      gridX=round(sqrt(maxexper));
      gridY=gridX;
      if gridX ~= sqrt(maxexper)
	  fprintf('WARNING: maxexp should be a square number!!!\n');
      end
      t=1:maxexper;
      [u,v]=ind2sub([gridX, gridY],t);
      pmultu=uval* (u - min(u))/(max(u)-min(u));
      pmultv=uval* (v - min(v))/(max(v)-min(v));
      params(1)=pmultu(jobid)*params(1);
      params(2)=pmultv(jobid)*params(2);
  else
      error('not more than 2 params at a time, please');
  end

  for k=1:length(keys)
      key=char(keys{k});
      %params = setfield(params,key,ini.GetValues(sec,key));
      opt = setfield(opt, key, params(k));
  end

  for k=1:length(keys)
      key=char(keys{k});
      %params = setfield(params,key,ini.GetValues(sec,key));
      opt = setfield(opt, key, params(k));
  end

  % write out new opt file
  status = writeDCOptions(opt,conffile);


  rng(1);

  allscensfile=fullfile(confdir,'doscens.txt');
  if ~exist(allscensfile)
	  warning('Doing the standard PETS TUD combo...');
	  dlmwrite(fullfile(confdir,'doscens.txt'),[23 25 27 71 72 42]);
  end
  allscen=dlmread(fullfile(confdir,'doscens.txt'));
  allscen

  learniter=str2double(learniter)

  mets2d=zeros(max(allscen),14);
  mets3d=zeros(max(allscen),14);  
  ens=zeros(max(allscen),5);

  for scenario=allscen
	  fprintf('jobid: %d,   learn iteration %d\n',jobid,learniter);
	  scensolfile=sprintf('%s/prt_res_%03d-scen%02d.mat',resdir,jobid,scenario)

	  try
	    load(scensolfile);
	  catch err
	    fprintf('Could not load result: %s\n',err.message);
	    [metrics2d, metrics3d, energies, stateInfo]=swDCTracker(scenario,conffile);
	    save(scensolfile,'stateInfo','metrics2d','metrics3d','energies');
	  end	  

	  mets2d(scenario,:)=metrics2d;
	  mets3d(scenario,:)=metrics3d;
	  ens(scenario,:)=double(energies);
	  infos(scenario).stateInfo=stateInfo;

  end


    
  save(resultsfile,'opt','mets2d','mets3d','ens','infos','allscen');
  
  
  % remove temp scene files
  for scenario=allscen
    scensolfile=sprintf('%s/prt_res_%03d-scen%02d.mat',resdir,jobid,scenario)
    if exist(scensolfile,'file')
      delete(scensolfile);
    end
  end
end



fprintf('Job done (%.2f min = %.2fh = %.2f sec per sequence)\n', ...
    toc(trainStartTime)/60,toc(trainStartTime)/3600,toc(trainStartTime)/numel(allscen));


% evaluate what we have so far
bestexper=combineResultsRemote(settingsDir,maxexper);

querystring=sprintf('qstat -t | grep %s | wc -l',settingsDir);
[rs,rjobs] = system(querystring); rjobs=str2double(rjobs)-1; % subtract currently running

fprintf('%d other jobs still running\n',rjobs);
