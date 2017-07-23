function mhsret=shrinkSplines(allpoints, mhs, used,labeling,T,E)
% create new trajectory hypotheses by srinking existing ones
% cf. Fig 10(b), PAMI


mhsnew=getEmptyModelStruct;

global opt
nMaxAddExt=opt.nMaxAddExt;
shrtail=4; shrhead=4;
nCurModels=length(mhs);
addedShrink=0;
nUsed=length(used);
% if nUsed>1
    for m=randperm(nUsed)
        if addedShrink>=nMaxAddExt, break; end
                
        mod1=used(m);
        supp=find(labeling==mod1);
        t=allpoints.tp(supp);
        [sortedsupp, sidx]=sort(t);
        if length(unique(sortedsupp))>=4+mean([shrtail shrhead])
            
            tnotail=t(sidx(1+shrtail:end));tr=tnotail;
            suppnotail=supp(sidx(1+shrtail:end));
                
%             tr=tnotail+opt.randFit*rand(1,length(tnotail));
            xy=[allpoints.xp(suppnotail);allpoints.yp(suppnotail)];
                    
            sfit=splinefit(tr,xy,mhs(used(m)).pieces);
            sfit=adjustSplineStruct(sfit, min(tnotail), max(tnotail), allpoints, T, 0, [], [], []);        
            addedShrink=addedShrink+1;
            
            mhsnew(addedShrink) = sfit;                   
            
            
            tnohead=t(sidx(1:end-shrhead));tr=tnohead;
            suppnohead=supp(sidx(1:end-shrhead));
                
%             tr=tnohead+opt.randFit*rand(1,length(tnohead));
            xy=[allpoints.xp(suppnohead);allpoints.yp(suppnohead)];
                    
            sfit=splinefit(tr,xy,mhs(used(m)).pieces);
            sfit=adjustSplineStruct(sfit, min(tnohead), max(tnohead), allpoints, T, 0, [], [], []);  
            addedShrink=addedShrink+1;
            mhsnew(addedShrink) = sfit;                   
            
                    
        end
                
    end
            
    
% end

allnewgoodness=[mhsnew.labelCost];
mhsnew=mhsnew(E>allnewgoodness);
for m=1:length(mhsnew), mhsnew(m).lastused=0; end

mhsret=mhs;
if ~isempty(mhsnew), mhsret=[mhs mhsnew]; end
end