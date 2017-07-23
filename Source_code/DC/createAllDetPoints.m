function alldpoints=createAllDetPoints(detections)
% put all data points (detections) into one struct
% 


T=size(detections,2); 
alldpoints.xp=[];alldpoints.yp=[];alldpoints.sp=[];alldpoints.tp=[];
if isfield(detections(1),'dirxi'), alldpoints.dirxi=[]; alldpoints.diryi=[]; end
if isfield(detections(1),'dirxw'), alldpoints.dirxw=[]; alldpoints.diryw=[]; end
if isfield(detections(1),'dirx'), alldpoints.dirx=[]; alldpoints.diry=[]; end
for t=1:T
    alldpoints.xp=[alldpoints.xp detections(t).xp];
    alldpoints.yp=[alldpoints.yp detections(t).yp];
    alldpoints.sp=[alldpoints.sp detections(t).sc];
    alldpoints.tp=[alldpoints.tp t*ones(1,length(detections(t).xp))];    
    if isfield(detections(t),'dirxi')
        alldpoints.dirxi=[alldpoints.dirxi detections(t).dirxi];
        alldpoints.diryi=[alldpoints.diryi detections(t).diryi];
    end
    if isfield(detections(t),'dirxw')
        alldpoints.dirxw=[alldpoints.dirxw detections(t).dirxw];
        alldpoints.diryw=[alldpoints.diryw detections(t).diryw];
    end
    if isfield(detections(t),'dirx')
        alldpoints.dirx=[alldpoints.dirx detections(t).dirx];
        alldpoints.diry=[alldpoints.diry detections(t).diry];
    end
end