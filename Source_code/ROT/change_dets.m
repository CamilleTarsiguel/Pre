for fr = 1:1194
    ind = M(:,1) == fr;
    detections(fr).x = M(ind,3);
    detections(fr).y = M(ind,4);
    detections(fr).w = M(ind,5);
    detections(fr).h = M(ind,6);
end
    