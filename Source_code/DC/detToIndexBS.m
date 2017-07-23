function [allt, allos, allos2, allkos]=detToIndexBS(splines,newLabeling, alldpoints)
%% find out, which point belongs to which spline piece (B-spline)

N=length(splines);
Npts=length(newLabeling);

splos=0;

allt=zeros(1,Npts);
allos=zeros(1,Npts);
allos2=zeros(1,Npts);
allkos=zeros(1,Npts);


knoffset=0;
for id=1:N
    
    % get spline
    spl=splines(id);
    
    % find points that belong to this spline
    idPts=find(newLabeling==id);       
    alldpts.tp=alldpoints.tp(idPts);
    
    % how many are there?
    npts=length(alldpts.tp);
        
%     breaks=spl.breaks;
    
    tr= alldpts.tp;
%     [~, index] = histc(tr,[-inf,spl.breaks(2:spl.pieces),inf]);
    [~, index] = histc(tr,[-inf,spl.bspline.knots(spl.bspline.order+1:spl.bspline.number),inf]);
%     splines(id).ptindex=index;
%     index=spl.ptindex;    
%     tr
%     index
%     breaks(index)
%     locc=tr-breaks(index);
    
    lcnt=0;
    
    sn=spl.bspline.number;
    
    for det=1:npts  
        lcnt=lcnt+1;
%         t=locc(det);       
        ilcnt=index(lcnt);
        bslcnt=ilcnt+3;
        
%         os = (index(lcnt)-1)*8 + splos;
        os=splos;
        
        allt(idPts(lcnt))=tr(det);
        allos(idPts(lcnt))=os+ilcnt-1;
        allos2(idPts(lcnt))=os+ilcnt-1+sn;
        kos=bslcnt-2 + knoffset;
        
        allkos(idPts(lcnt))=kos-1;
        
        
        
    end
    splos = splos + sn*2; 
    knoffset = knoffset + length(spl.bspline.knots);
end