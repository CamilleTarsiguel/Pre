function [allt, allos]=detToIndex(splines,newLabeling, alldpoints)
%% find out, which point belongs to which spline piece (PP)

N=length(splines);
Npts=length(newLabeling);

splos=0;

allt=zeros(1,Npts);
allos=zeros(1,Npts);

for id=1:N
    
    % get spline
    spl=splines(id);
    
    % find points that belong to this spline
    idPts=find(newLabeling==id);       
    alldpts.tp=alldpoints.tp(idPts);
    
    % how many are there?
    npts=length(alldpts.tp);
        
    breaks=spl.breaks;
    
    tr= alldpts.tp;
    [~, index] = histc(tr,[-inf,spl.breaks(2:spl.pieces),inf]);
%     splines(id).ptindex=index;
%     index=spl.ptindex;    
%     tr
%     index
%     breaks(index)
    locc=tr-breaks(index);
    
    lcnt=0;
    
    for det=1:npts  
        lcnt=lcnt+1;
        t=locc(det);       
        
        os = (index(lcnt)-1)*8 + splos;
        allt(idPts(lcnt))=t;
        allos(idPts(lcnt))=os;
        
        
        
    end
    splos = splos + spl.pieces*8; 
end