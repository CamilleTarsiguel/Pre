function [mhsret splitPairs pts]=splitSplines(alldpoints, splines, used,labeling,T,E)
% create new trajectory hypotheses by splitting existing ones
% along no-detection gaps, Fig. 10(c), PAMI
%

% TODO

mhsnew=getEmptyModelStruct();

global opt sceneInfo
tau=sceneInfo.targetSize/2;

nMaxAddExt=opt.nMaxAddExt; % TODO change parameter
nCurModels=length(splines);
addedExt=0;
nUsed=length(used);
stateInfo=getStateFromSplines(splines, struct('F',T));

splitPairs=[];
pts=[];
for m=used
    mhs=splines(m);
    tt=mhs.start:mhs.end;
    splinexy=[stateInfo.X(:,m)';stateInfo.Y(:,m)'];
    
    
    detsWithin=-1*ones(1,T); % detections close to model
    insideAL=zeros(1,T); % inside areaLimits
    
    %% data fidelity
    for t=tt
        sx=splinexy(1,t); sy=splinexy(2,t);
        
        detind=find(alldpoints.tp==t);
        insideAL(t)=1;
        detsWithin(t)=0;
        
        sp=[sx*ones(1,numel(detind)); sy*ones(1,numel(detind))];
        dp=[alldpoints.xp(detind);alldpoints.yp(detind)];
        d=sp-dp;
        
        alldist=sqrt(sum(d.^2));
        
        inthisframe=alldist < tau;
        
        detsWithin(t)=sum(inthisframe);
        
    end
    nodetframes=[0 ~detsWithin(insideAL>0) 0];
    c=find(diff(nodetframes)==1); d=find(diff(nodetframes)==-1);
    
    occgaps=(d-c);
    
    
    %     nodetframes
    [maxgap mgpos]=max(occgaps);
    %
    %     maxgap
    %     c
    %     d
    %     mgpos
    safmar=3;
    if ~isempty(maxgap)
        if maxgap>=5 && c(mgpos)>safmar && d(mgpos)+safmar<length(tt)
            tt1=tt(1):tt(c(mgpos)-safmar);
            tt2=tt(d(mgpos)+safmar):tt(end);
            
            
            
            %     tt1
            %     tt2
            
            drawsfit=[];
            splshr=0;
            % fit through splines
            if length(tt1)>=4
                xy=ppval(mhs,tt1); %pts on spline
                tr=tt1;
%                 tr=tt1+opt.randFit*rand(1,length(tt1));
                addedExt=addedExt+1;
                
                sfit=splinefit(tr,xy,round(max(opt.minCPs,(tr(end)-tr(1))*opt.ncpsPerFrame )));
                sfit=adjustSplineStruct(sfit, min(tt1), max(tt1), alldpoints, T, 0, [], [], []); 
                mhsnew(addedExt) = sfit;
                drawsfit=[drawsfit sfit];
                splshr=splshr+1;
            end
            
            if length(tt2)>=4
                xy=ppval(mhs,tt2); %pts on spline
                tr=tt2;
%                 tr=tt2+opt.randFit*rand(1,length(tt2));
                addedExt=addedExt+1;
                
                sfit=splinefit(tr,xy,round(max(opt.minCPs,(tr(end)-tr(1))*opt.ncpsPerFrame )));
                sfit=adjustSplineStruct(sfit, min(tt2), max(tt2), alldpoints, T, 0, [], [], []); 
                mhsnew(addedExt) = sfit;
                drawsfit=[drawsfit sfit];
                splshr=splshr+1;
                
            end
%             drawsfit=[mhs drawsfit];
%             length(drawsfit)
%             splshr
%             prepFigure; drawSplines(drawsfit,1:length(drawsfit),0,[],1:T)
%             pause
            
            
            if splshr==2
                splitPairs(end+1,:)=[m addedExt-1 addedExt];
                pt1=find(alldpoints(:).tp>=tt1(1) & alldpoints(:).tp<=tt1(end));
                pt2=find(alldpoints(:).tp>=tt2(1) & alldpoints(:).tp<=tt2(end));
                lab=find(labeling==m);  
                pts(end+1).pt1=intersect(pt1,lab);
                pts(end).pt2=intersect(pt2,lab);
%                 [gd sgc0]=getSplineGoodness(mhs,1,alldpoints,T);
%                 sgc0
%                 sgc1
%                 sgc2
%                 sum(sgc0)
%                 sum([sgc1 sgc2])

            end
            
        end
        
        
    end
    %     d=[d length(nodetframes)]+1;
    
    %     for o=reloccgaps
    %         c(o):c(o+1)
    %         1:c(o)
    %         d(o)+1:length(nodetframes)-1)
    %     end
    %     pause
end

allnewgoodness=[mhsnew.labelCost];
mhsnew=mhsnew(E>allnewgoodness);
for m=1:length(mhsnew), mhsnew(m).lastused=0; end


mhsret=splines;
% length(mhsret)
if ~isempty(mhsnew), mhsret=[splines mhsnew]; end

% length(mhsret)
end