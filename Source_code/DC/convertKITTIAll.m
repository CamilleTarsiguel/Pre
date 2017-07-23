%%
for s=0:20
    cd /home/amilan/storage/databases/KITTI/tracking/devkit_tracking/matlab
    tracklets=readLabels('../../training/label_02',s);
    cd /home/amilan/research/projects/dctracking
    convertKITTIToCVML(tracklets,sprintf('/home/amilan/storage/databases/KITTI/tracking/training/label_02/%04d.xml',s),{'Car','Van'});
end