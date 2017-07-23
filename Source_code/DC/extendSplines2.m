function mhsret=extendSplines2(alldpoints, mhs, used,opt)
% create new trajectory hypotheses by extending existing ones
%
% all used ones are extended 2 frames back and forward


mhsnew=getEmptyModelStruct;

exttail=2;
exthead=2;
nCurModels=length(mhs);
addedExt=0;
nUsed=length(used);

% length(used)
for m=used
    splnew=mhs(m);
%     splnew
    if splnew.start>exttail
        splnew.start=splnew.start-exttail;
        addedExt=addedExt+1;
        
        [splnew.labelCost, splnew.lcComponents] = ...
            getSplineGoodness(splnew,1,alldpoints,opt.seqLength);
        mhsnew(addedExt) = splnew;
%         fprintf('added 1\n');
    end
    
    if splnew.end <= opt.seqLength-exthead
        splnew=mhs(m);
        
        splnew.end=splnew.end+exthead;
        addedExt=addedExt+1;
        
        [splnew.labelCost, splnew.lcComponents] = ...
            getSplineGoodness(splnew,1,alldpoints,opt.seqLength);
        mhsnew(addedExt) = splnew;
%         fprintf('added 2\n');
    end
    
    if splnew.start>exttail && splnew.end <= opt.seqLength-exthead
        splnew=mhs(m);
        splnew.end=splnew.end+exthead;
        splnew.start=splnew.start-exttail;
        addedExt=addedExt+1;
        
        [splnew.labelCost, splnew.lcComponents] = ...
            getSplineGoodness(splnew,1,alldpoints,opt.seqLength);
        mhsnew(addedExt) = splnew;
%         fprintf('added 3\n');
    end
end

% length(mhsnew)
% pause
for m=1:length(mhsnew), mhsnew(m).lastused=0; end

mhsret=mhs;
if ~isempty(mhsnew), mhsret=[mhs mhsnew]; end
end