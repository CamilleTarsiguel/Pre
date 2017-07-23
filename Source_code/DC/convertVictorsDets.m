allscen=[1,3,4,6,7,8,10,12];
allscen=[0,2,5,11,13,15,17,19];
minx=Inf; maxx=-Inf; miny=Inf; maxy=-Inf;
for s=allscen
    %load file
    load(sprintf('/home/amilan/research/projects/dctracking/tmp/TestingDetections/S%02d/detections.mat',s));

    detections=[];
    gtInfo=[];
    % for each frame
    for t=1:length(detections_gt)
        
        
        
        
        % for each det in frame
        for n=1:size(detections_gt{t}.bbox,1)
            bb=detections_gt{t}.bbox(n,1:4); % center x,y, width, height
            ct=detections_gt{t}.centroids(n,:) * 1000; % x,y plane and height in world mm
            xp=ct(1); yp=ct(2);
            
            if xp<minx, minx=xp; end
            if xp>maxx, maxx=xp; end
            if yp<miny, miny=yp; end
            if yp>maxy, maxy=yp; end
            
            w=bb(3); h=bb(4);
            x1=bb(1)-w/2; y1=bb(2)-h/2;
            
            detections(t).bx(n)=x1;
            detections(t).by(n)=y1;
            detections(t).xp(n)=bb(1);
            detections(t).yp(n)=bb(2)+h/2;
            detections(t).xw(n)=ct(1);
            detections(t).yw(n)=ct(2);            
            detections(t).ht(n)=h;
            detections(t).wd(n)=w;
            detections(t).sc(n)=1;
            detections(t).xi(n)=bb(1);
            detections(t).yi(n)=bb(2)+h/2;
            
            id=detections_gt{t}.GTAssoc(n)+1;
            gtInfo.Xgp(t,id)=ct(1);gtInfo.Ygp(t,id)=ct(2);
            gtInfo.W(t,id)=w;gtInfo.H(t,id)=h;
            gtInfo.Xi(t,id)=bb(1);gtInfo.Yi(t,id)=bb(2)+h/2;
        end
    end
    t=t+1;
    
    % Victor is missing last frame
    detections(t).bx = [];
    
    gtInfo.Xgp(t,:)=0;gtInfo.Ygp(t,:)=0;
    gtInfo.W(t,:)=0;gtInfo.H(t,:)=0;
    gtInfo.Xi(t,:)=0;gtInfo.Yi(t,:)=0;
%     gtInfo.X=gtInfo.Xgp;gtInfo.Y=gtInfo.Ygp; % 3d
    gtInfo.X=gtInfo.Xi;gtInfo.Y=gtInfo.Yi; % 3d
    
    gtInfo.frameNums=0:t-1;
    
    save(sprintf('/home/amilan/storage/databases/KITTI/tracking/training/det_02/Victor/%04d.mat',s),'detections');
    save(sprintf('/home/amilan/storage/databases/KITTI/tracking/training/label_02/Victor/%04d.mat',s),'gtInfo');
end