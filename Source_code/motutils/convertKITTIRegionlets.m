%%  convertKITTIRegionlets

detfolder='/home/amilan/storage/databases/KITTI/tracking/testing/det_02/regionlets';

classes={'Pedestrian'};
detstr='-peds';
detthr=0;
seq=0;
seqdir=fullfile(detfolder,sprintf('%04d',seq));
while exist(seqdir,'dir')
    
    fr=0;
    framefile=fullfile(seqdir,sprintf('%06d.txt',fr));
    
    
    
    clear C detections
    while exist(framefile,'file')
        
        t=fr+1;
        fid = fopen(framefile);
        
        C = textscan(fid, '%s %d %d %d %f %f %f %f %d %d %d %d %d %d %d %f');
        nlines=size(C{1},1);
        
        bx=[];by=[];bw=[];bh=[];xi=[];yi=[];sc=[];
        
        for d=1:nlines
            %detection too weak
            if C{16}(d) < detthr, continue; end
            
            % class ignored
            if ~ismember(C{1}(d),classes), continue; end
            
            
            bbox=[C{5}(d) C{6}(d) C{7}(d) C{8}(d) C{16}(d)] + 1; % +1 because KITTI 0-based
            bw_=bbox(3)-bbox(1);bh_=bbox(4)-bbox(2);
            sc_=1/(exp(-bbox(5))+1);
            
            bx = [bx bbox(1)];            by = [by bbox(2)];
            bw = [bw  bw_];            bh = [bh bh_];
%             xi = [xi bbox(1)+bw_/2];
%             yi = [yi bbox(2)+bh_];
            sc = [sc sc_];
            
        end
        
        detections(t).bx = bx;
        detections(t).by = by;
        detections(t).wd = bw;
        detections(t).ht = bh;
        detections(t).xi = bx+bw./2;
        detections(t).yi = by+bh;
        detections(t).sc = sc;
        
        detections(t).xp=detections(t).xi;
        detections(t).yp=detections(t).yi;
        
%         
        fclose(fid);
        
        fr=fr+1;
        framefile=fullfile(seqdir,sprintf('%06d.txt',fr));
%         break
    end
    
    % save det file
    detfile=fullfile(detfolder,sprintf('%04d%s',seq,detstr));
    save(detfile,'detections');
%     break
    seq=seq+1;
    seqdir=fullfile(detfolder,sprintf('%04d',seq));
end

