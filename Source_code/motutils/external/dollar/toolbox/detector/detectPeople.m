function detFile=detectPeople(sceneInfo)
% run P. Dollar's pedestrian detector.




imgFolder = sceneInfo.imgFolder;

% detection's folder and file
detFile = fullfile(sceneInfo.imgFolder,'det.txt');
if exist(detFile,'file')
    fprintf('Detections found: %s\n',detFile);
    return;
end


load('models/AcfInriaDetector.mat')

detScale=1;
% detScale=0.6;

blowUp=1;
fprintf('RESCALING DETECTOR: %f\n',detScale);
detector = acfModify(detector,'rescale',detScale);


imgExt = sceneInfo.imgExt;

imgMask=[imgFolder,'/*',imgExt];
dirImages = dir(imgMask);


F=length(dirImages);
%         F=10;
filecells=cell(1,F);
for t=1:F
    filecells{t} = [imgFolder,dirImages(t).name];
end

fprintf('Detecting %s (%d frames)\n',sceneInfo.sequence,F);


delete(detFile);

% detect all images
%         detector.opts.pNms.type='none';
bbx = acfDetect(filecells,detector);

% write out
writeDets(bbx,detFile);

end
