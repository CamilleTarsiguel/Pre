%% Script pour ACF Modify


pNms = struct('type', 'max', 'overlap', 0.3, 'ovrDnm', 'min');
pModify=struct('pNms', pNms3);
detector6 = acfModify(detector, pModify);
BB6 = acfDetect('~/Documents/MATLAB/Pre/data/SMSQ20/frame-0.jpg', detector6);
bb = BB6(BB6(:,5)>100,:)
I6 = insertObjectAnnotation(I, 'rectangle', bb(:,1:4),bb(:,5));
imshow(I6);
