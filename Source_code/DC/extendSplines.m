function mhsret=extendSplines(allpoints, mhs, used,labeling,T,E)
% create new trajectory hypotheses by extending existing ones
% 


mhsnew=getEmptyModelStruct;

global opt
nMaxAddExt=opt.nMaxAddExt;
exttail=2;
exthead=2;
nCurModels=length(mhs);
addedExt=0;
nUsed=length(used);

% if nUsed>1
    for m=randperm(nUsed)
        if addedExt>=nMaxAddExt, break; end
        
        % extend linearly
        mod1=used(m);
        supp=find(labeling==mod1);
        t=allpoints.tp(supp);
        [sortedsupp sidx]=sort(t);
        if length(unique(sortedsupp))>=4
            
            xy=[allpoints.xp(supp);allpoints.yp(supp)];
            
            extended=0;
            if sortedsupp(1)-exttail>0
                
                ntail=sidx(1:4);
                tailpts=supp(ntail);
                xytail=[allpoints.xp(tailpts);allpoints.yp(tailpts)];
                ttail=allpoints.tp(tailpts);
                if length(unique(ttail))>2
                    tailline=splinefit(ttail,xytail,1,2);
                    taillinepts=ppval(tailline,sortedsupp(1)-exttail:sortedsupp(1)-1);
                    
                    xy=[taillinepts xy];
                    t=[sortedsupp(1)-exttail:sortedsupp(1)-1 t];tr=t;
%                     tr=t+opt.randFit*rand(1,length(t));
                    
                    
                    sfit=splinefit(tr,xy,mhs(used(m)).pieces);                                        
                    sfit=adjustSplineStruct(sfit, min(t), max(t), allpoints, T, 0, [], [], []);
                    addedExt=addedExt+1;
                    mhsnew(addedExt) = sfit;
                    
                    extended=1;
                end
                
            end
            
            if sortedsupp(end)+exthead<T
                t=allpoints.tp(supp);
                xy=[allpoints.xp(supp);allpoints.yp(supp)];
                
                nhead=sidx(end-3:end);
                headpts=supp(nhead);
                xyhead=[allpoints.xp(headpts);allpoints.yp(headpts)];
                thead=allpoints.tp(headpts);
                if length(unique(thead))>2
                    headline=splinefit(thead,xyhead,1,2);
                    headlinepts=ppval(headline,sortedsupp(end)+1:sortedsupp(end)+exthead);
                    
                    xy=[xy headlinepts];
                    t=[t sortedsupp(end)+1:sortedsupp(end)+exthead];tr=t;
%                     tr=t+opt.randFit*rand(1,length(t));
                    
                    
                    sfit=splinefit(tr,xy,mhs(used(m)).pieces);
                    sfit=adjustSplineStruct(sfit, min(t), max(t), allpoints, T, 0, [], [], []);                                                            
                    addedExt=addedExt+1;
                    mhsnew(addedExt) = sfit;
                    
                    extended=1;
                end
            end
%             
%             if extended
%             end
        end
        
      
%         mhsnew
    end
    
% end


allnewgoodness=[mhsnew.labelCost];
mhsnew=mhsnew(E>allnewgoodness);
for m=1:length(mhsnew), mhsnew(m).lastused=0; end

mhsret=mhs;
if ~isempty(mhsnew), mhsret=[mhs mhsnew]; end
end