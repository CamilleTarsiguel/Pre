function sceneInfo=getSceneInfo(scenario, opt)
% fill all necessary information about the
% scene into the sceneInfo struct
%
% Required:
%   detfile         detections file (.idl or .xml)
%   frameNums       frame numbers (eg. frameNums=1:107)
%   frameRate       frame rate of the video sequence (default: 25)
%   vidFile         video file
%   targetSize      approx. size of targets (default: 5 on image, 350 in 3d)
%
% Required for 3D Tracking only
%   trackingArea    tracking area
%   camFile         camera calibration file (.xml PETS format)
%
% Optional:
%   gtFile          file with ground truth bounding boxes (.xml CVML)
%   initSolFile     initial solution (.xml or .mat)
%   targetAR        aspect ratio of targets on image
%   bgMask          mask to bleach out the background
%   dataset         name of the dataset
%   sequence        name of the sequence
%   scenario        sequence number



if nargin<2
    global opt
end

dbfolder='data';
sceneInfo.dataset='TUD';
sceneInfo.sequence='TUD-Stadtmitte';
sceneInfo.frameNums=7022:7200;
sceneInfo.frameRate=25;

sceneInfo.imgFolder=fullfile(dbfolder,filesep,'cvpr10_tud_stadtmitte',filesep);
sceneInfo.imgFileFormat='DaMultiview-seq%04d.png';
% image dimensions
[sceneInfo.imgHeight, sceneInfo.imgWidth, ~]= ...
    size(imread([sceneInfo.imgFolder sprintf(sceneInfo.imgFileFormat,sceneInfo.frameNums(1))]));

% detections
sceneInfo.detfile=fullfile(dbfolder,'det',sceneInfo.dataset,'TUD-stadtmitte-det.xml');

%% tracking area
% if we are tracking on the ground plane
% we need to explicitly secify the tracking area
% otherwise image = tracking area
sceneInfo.trackingArea=[-19, 12939, -48, 10053];

% camera calibration
cameraconffile=fullfile(dbfolder,filesep,'TUD-Stadtmitte-calib.xml');
sceneInfo.camFile=cameraconffile;

if ~isempty(sceneInfo.camFile)
    sceneInfo.camPar=parseCameraParameters(sceneInfo.camFile);
    %%
end

%% target size
sceneInfo.targetSize=20;                % target 'radius'
sceneInfo.targetSize=sceneInfo.imgWidth/30;
if opt.track3d, sceneInfo.targetSize=350; end

%% target aspect ratio
sceneInfo.targetAR=1/3;


%% ground truth
sceneInfo.gtFile=fullfile(dbfolder,'gt',sceneInfo.dataset,[sceneInfo.sequence '.xml']);
fprintf('GT File: %s\n',sceneInfo.gtFile);

global gtInfo
sceneInfo.gtAvailable=0;
if ~isempty(sceneInfo.gtFile)
    sceneInfo.gtAvailable=1;
    % first determine the type
    [pathtogt, gtfile, fileext]=fileparts(sceneInfo.gtFile);
    
    if strcmpi(fileext,'.xml') % CVML
        gtInfo=parseGT(sceneInfo.gtFile);
    elseif strcmpi(fileext,'.mat')
        % check for the var gtInfo
        fileInfo=who('-file',sceneInfo.gtFile);
        varExists=0; cnt=0;
        while ~varExists && cnt<length(fileInfo)
            cnt=cnt+1;
            varExists=strcmp(fileInfo(cnt),'gtInfo');
        end
        
        if varExists
            load(sceneInfo.gtFile,'gtInfo');
        else
            warning('specified file does not contain correct ground truth');
            sceneInfo.gtAvailable=0;
        end
    end
    
    if opt.track3d
        if ~isfield(gtInfo,'Xgp') || ~isfield(gtInfo,'Ygp')
            [gtInfo.Xgp, gtInfo.Ygp]=projectToGroundPlane(gtInfo.X, gtInfo.Y, sceneInfo);
        end
    end
    gtInfo.Xi=gtInfo.X; gtInfo.Yi=gtInfo.Y;
    
    if opt.remOcc
        gtInfo=removeOccluded(gtInfo);
    end
    
        if strcmpi(fileext,'.xml'),     save(fullfile(pathtogt,[gtfile '.mat']),'gtInfo'); end
end

%% cut GT to tracking area
if  sceneInfo.gtAvailable && opt.track3d && opt.cutToTA
    gtInfo=cutGTToTrackingArea(gtInfo,sceneInfo);
end
%% check
if opt.track3d
    if ~isfield(sceneInfo,'trackingArea')
        error('tracking area [minx maxx miny maxy] required for 3d tracking');
    elseif ~isfield(sceneInfo,'camFile')
        error('camera parameters required for 3d tracking. Provide camera calibration or siwtch to 2d tracking');
    end
end

%% shift target center from foot position to center of BB?
sceneInfo.yshift=0;

sceneInfo.scenario=scenario;
end
