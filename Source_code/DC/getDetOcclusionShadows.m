function occmap=getDetOcclusionShadows(tt)
% compute a mask for detections overlays
% (not used)

global detections sceneInfo
occfilename=[sceneInfo.imgFolder sprintf('occmaps-%i-%i.mat',sceneInfo.frameNums(1),sceneInfo.frameNums(end))];
if exist(occfilename,'file')
    load(occfilename)
    return;
end

F=1;
if nargin<1
    F=length(detections);
    tt=1:F;
end


% occmap=false(sceneInfo.imgHeight,sceneInfo.imgWidth,F);

for t=tt
%     clf
%     t
%     detections(t)
    occmap(t).map=false(sceneInfo.imgHeight,sceneInfo.imgWidth);
    ndets=length(detections(t).xp);
    for d=1:ndets
        xleft=round(detections(t).bx(d));
        ytop =round(detections(t).by(d));
        xright=round(xleft+detections(t).wd(d));
        ybot=round(ytop+detections(t).ht(d));
        [xleft xright ytop ybot]= ...
    clampBBox(xleft, xright, ytop, ybot, sceneInfo.imgWidth, sceneInfo.imgHeight);
        
        
%         [ytop ybot xleft xright]
%         pause
        occmap(t).map(ytop:ybot,xleft:xright)=1;
%         occmap(ytop:ybot,xleft:xright,t)=1;
        % scene occluder
        occmap(t).map(41:400,405:455)=1;
    end
%     occtoshow=occmap(:,:,t);
%     size(occtoshow)
%     imshow(occmap(:,:,t));
%     pause
    
end

save(occfilename,'occmap');
end