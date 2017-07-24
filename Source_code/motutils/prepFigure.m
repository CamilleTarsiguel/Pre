function prepFigure()
% prepare figure for showing state
%

global sceneInfo opt scenario detections;

% figh=findobj('type','figure','name','optimization');

% if isempty(figh), figh=figure('name','optimization'); end
% set(figh);

clf;
hold on;
box on

if ~opt.track3d
    set(gca,'Ydir','reverse');
end

% if trackingArea not defined, used min max from detections
if isfield(sceneInfo,'trackingArea')
    trAr=sceneInfo.trackingArea;
else
    trAr=[ min([detections(:).xp]), ...
        max([detections(:).xp]), ...
        min([detections(:).yp]), ...
        max([detections(:).yp])];
end

xlim(trAr(1:2))
ylim(trAr(3:4))
if ~opt.track3d
    ylim([sceneInfo.imTopLimit trAr(4)]);
end

zlim([0 length(sceneInfo.frameNums)+1])

view(-78,4)
if ~opt.track3d
    view(8,12);
end


F=length(sceneInfo.frameNums);
for t=1:F
    ndet=length(detections(t).bx);
    for d=1:ndet
        x1=detections(t).bx(d); y1=detections(t).by(d);
        y2=y1+detections(t).ht(d); x2=x1+detections(t).wd(d);
        sc=detections(t).sc(d);
        
        %     line([x1 x1 x2 x2 x1],[y1 y2 y2 y1 y1],[t t t t t],'linewidth',sc*5);
    end
end
end