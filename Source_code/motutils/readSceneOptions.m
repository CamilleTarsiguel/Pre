function sceneInfo=readSceneOptions(inifile)
% parse sceneInfo

if ~exist(inifile,'file')
    fprintf('WARNING! Scene file %s does not exist! Using TUD Campus ...\n',inifile);
    inifile='scenes/scene2D.ini';
end


ini=IniConfig();
ini.ReadFile(inifile);

sceneInfo=[];

% Sequence Information
sceneInfo=fillInOpt(sceneInfo,ini,'Source');

% Additional Info
%sceneInfo=fillInOpt(sceneInfo,ini,'Miscellaneous');

sceneInfo.scenario=0;

sceneInfo=checkScene(sceneInfo);

end


function sceneInfo = fillInOpt(sceneInfo, ini, sec)
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
    sceneInfo = setfield(sceneInfo,key,val);
end

end

function sceneInfo=checkScene(sceneInfo)
% check sceneInfo for correctness
%
% TODO: Insert on-the-fly detector

% mandatory fields
requiredFields={'imgFolder','frameRate'};

for f=requiredFields
    fchar=char(f);
    assert(isfield(sceneInfo,fchar), ...
        ['Field ' fchar ' is required']);
end
reader = VideoReader(sceneInfo.imgFolder);
% sceneInfo.imgFolder=[sceneInfo.imgFolder,'/'];
% 
% assert(exist(sceneInfo.imgFolder,'dir')>0, ...
%     ['Image Folder ' sceneInfo.imgFolder ' does not exist.']);

% figure out frames, format, etc.
sceneInfo.frameNums = 1:reader.numberOfFrames ;

% append file extension
%sceneInfo.imgFileFormat=[sceneInfo.imgFileFormat, sceneInfo.imgExt];

% image dimensions
sceneInfo.imgHeight = reader.Height;
sceneInfo.imgWidth = reader.Width;

% run detector if no detfile given
if ~isfield(sceneInfo,'detfile')
    detfile=detectPeople(sceneInfo);
    sceneInfo.detfile = fullfile(detfile);

end

if ~isfield(sceneInfo,'dataset')
    sceneInfo.dataset='unknown';
end
if ~isfield(sceneInfo,'sequence')
    sceneInfo.sequence='unknown';
end

assert(exist(sceneInfo.detfile,'file')>0, ...
    ['Detection ' sceneInfo.detfile ' does not exist.']);



% ground truth
sceneInfo.gtAvailable=0;
if isfield(sceneInfo,'gtFile')
    if exist(sceneInfo.gtFile,'file')
        sceneInfo.gtAvailable=1;
        
        global gtInfo
        gtInfo=convertTXTToStruct(sceneInfo.gtFile);
        gtInfo.X=gtInfo.Xi; gtInfo.Y=gtInfo.Yi;
        
        
    end
end

%
if isfield(sceneInfo,'camFile')
    sceneInfo.camPar = parseCameraParameters(sceneInfo.camFile);
    if ~isfield(gtInfo,'Xgp') || ~isfield(gtInfo,'Ygp')
        [gtInfo.Xgp, gtInfo.Ygp]=projectToGroundPlane(gtInfo.X, gtInfo.Y, sceneInfo);
        gtInfo.X=gtInfo.Xgp;gtInfo.Y=gtInfo.Ygp;
    end       
end

end