function opt=readDCOptions(inifile)
% parse configuration for DCO Tracking

if ~exist(inifile,'file')
    fprintf('WARNING! Config file %s does not exist! Using default setting...\n',inifile);
    inifile='config/default2d.ini';
end


ini=IniConfig();
ini.ReadFile(inifile);

opt=[];

% General
opt=fillInOpt(opt,ini,'General');
% take care of frames vector
if isfield(opt,'ff') && isfield(opt,'lf')
    opt.frames=opt.ff:opt.lf;
    opt=rmfield(opt,'ff');opt=rmfield(opt,'lf');
end

% Hypothesis Space management
opt=fillInOpt(opt,ini,'Hypothesis Space');

% Initialization (EKF, DP, ...)
opt=fillInOpt(opt,ini,'Initialization');

% Detections
opt=fillInOpt(opt,ini,'Detections');
if isfield(opt,'sigA') && isfield(opt,'sigB')
    opt.detScale.sigA=opt.sigA;
    opt.detScale.sigB=opt.sigB;
    opt=rmfield(opt,'sigA');opt=rmfield(opt,'sigB');
end

% Functions
opt=fillInOpt(opt,ini,'Functions');

% Main Parameters
opt=fillInOpt(opt,ini,'Parameters');

% Optimization details
opt.disOpt=[];
opt.disOpt=fillInOpt(opt.disOpt,ini,'Discrete Optimization');
% opt.disOpt=disOpt;

% cont. optimization
opt.conOpt=[];
opt.conOpt=fillInOpt(opt.conOpt,ini,'Continuous Optimization');


opt=fillInOpt(opt,ini,'Miscellaneous');
end


function opt = fillInOpt(opt, ini, sec)
% loop through all keys in section and
% append to struct

keys = ini.GetKeys(sec);
for k=1:length(keys)
    key=char(keys{k});
    val=ini.GetValues(sec,key);
    
    % parameters are numeric
    if isstr(val) && strcmpi(sec,'Parameters')
        val=str2double(val);
    end
    opt = setfield(opt,key,val);
end

end
