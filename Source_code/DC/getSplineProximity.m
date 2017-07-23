function [prox, proxt, proxcost, merpairs]=getSplineProximity(mh,T,stInfo,opt)
% compute pairwise distances along all splines 
% Eq. (22), PAMI
%

F=length(mh);


% global opt sceneInfo
pF=opt.proxcostFactor;
proxcost=zeros(F);
prox=Inf*ones(F);
proxt=zeros(F);
merpairs=[];
if ~pF
    return;
end


% pairwise distance between splines

tt=1:T;

% compute all spline positions
if nargin<3,    stInfo=[]; end

if isempty(stInfo)
    X=zeros(T,F);    Y=zeros(T,F);
    for m1=1:F
        tt=mh(m1).start:mh(m1).end;
        splinexy=ppval(mh(m1),tt);
        X(tt,m1)=splinexy(1,:)';
        Y(tt,m1)=splinexy(2,:)';
    end
else
    X=stInfo.X; Y=stInfo.Y;
end
Xo=X; Yo=Y;

%
% isactive=false(1,F);
% for m=1:F
%     if any(labeling==m), isactive(m)=1;end
% end

% for m1=1:F
%     rounded(m1).breaks=round(mh(m1).breaks);
% end

% for m1=1:F
%     s1=mh(m1).start; e1=mh(m1).end;
%     s1=find(X(:,m1),1,'first');e1=find(X(:,m1),1,'last');
%     for m2=m1+1:F
%         s2=mh(m2).start; e2=mh(m2).end;
%         s2=find(X(:,m2),1,'first');e2=find(X(:,m2),1,'last');
%
%         mus=max(s1,s2);        mue=min(e1,e2);
%         if mus<=mue
%
%             evalPerFr=1;
%             tt=mus:mue;
%             %%%% CAREFUL! only for integer tt!!!
%             splinexy1=[X(tt,m1)';Y(tt,m1)'];
%             splinexy2=[X(tt,m2)';Y(tt,m2)'];
%
%             d=splinexy1-splinexy2;
%             squares=(d.^2);
%
%             sqrts=sqrt(sum(squares));
%             distalongtraj=sqrts/1;
%             [mindist atframe]=min(distalongtraj);
%
% %             [m1 m2]
% %             distalongtraj
%             prox(m1,m2)=mindist;
%             proxt(m1,m2)=mus+atframe/evalPerFr-1;
%         end
%     end
%
% end

% vectorized version
X(~X)=Inf; Y(~Y)=Inf;

prox2=Inf*ones(F,F,T);
for t=1:T
    tmp=repmat(X(t,:),F,1);
    ax=tmp(:)';
    tmp=tmp';
    bx=tmp(:)';
    
    tmp=repmat(Y(t,:),F,1);
    ay=tmp(:)';
    tmp=tmp';
    by=tmp(:)';
    
    dx=(ax-bx).^2;
    dy=(ay-by).^2;
    d=sqrt(dx+dy);
    
    % inf-inf = NaN
%     if ~isempty(find(isnan(d)))
%         d
%         ax
%         ay
%     end
    d(isnan(d))=Inf;
    
%     [size(prox2) F]
% if t==14
%     size(d)
%     t
%     fprintf('%.15f\n',ax);
%     fprintf('%.15f\n',bx);
%     fprintf('%.15f\n',ay);
%     fprintf('%.15f\n',by);
%     fprintf('%.15f\n',dx);
%     fprintf('%.15f\n',dy);
%     fprintf('%.15f\n',d);
% ax
% ay
% bx
% by
%     dx
%     dy
%      d
% end
%     pause
    prox2(:,:,t)=reshape(d,F,F);
end
% prox2
proxfromloop=prox2;
[prox2, proxt2]=min(prox2,[],3);
% if size(proxfromloop)>2
% [proxfromloop(1,2,:)]
% end
[prox3]=pF*sum(1./proxfromloop,3);
prox3=triu(prox3);

prox2=triu(prox2); prox2(isinf(prox2))=Inf;
prox2(~~(tril(ones(F))))=Inf;
% prox
% prox2
prox=prox2;
% assert(isequal(prox,prox2))
% assert(isequal(proxt,proxt2))

% prox(~proxt)=Inf;
% prox
% proxcost=sum(1./prox);
% proxcost=max(1./prox + 1./prox');
% proxcost=max(1./prox);
% prox(prox<200)=0;
% prox(prox<200)=1;
% prox(prox>=200)=Inf;

% sigi
siga=opt.conOpt.enParEexc(1);
sigb=(opt.conOpt.enParEexc(2)/2)*siga;
sigscale=pF;
proxsigi=sigscale-sigscale*1./(1+exp(-siga*proxfromloop+sigb));
% fprintf('%.15f\n',proxsigi(1,2,:)/100);
% proxfromloop
%%%% Tukey %%%%
% shift=sceneInfo.targetSize/4;
% k=sceneInfo.targetSize/2;
% ppp=proxfromloop(:)';
% 
% allxgrkind= abs(ppp-shift)>k;
% allxle0ind= ppp<shift;
% 
% y=(k^2/6)*(1-(1-((ppp-shift)./k).^2).^3);
% y(allxgrkind)=k^2/6;
% y(allxle0ind)=0;
% y=6/k^2*y;
% y=1-y;
% y=pF*y;
% 
% proxsigi=reshape(y,size(proxfromloop));

% sum(proxsigi(:))
% sum(proxsigi2(:))
% pause
proxsigimax=max(proxsigi,[],3);
proxsigisum=sum(proxsigi,3);
% proxsigisum
proxsigimax(1:F+1:F*F)=0;
proxsigisum(1:F+1:F*F)=0;
% proxsigisum
proxsigimax=triu(proxsigimax);
proxsigisum=triu(proxsigisum);
% proxsigisum

% proxsigi
% pause

proxcost=pF*1./prox;

prox3(1:F+1:F*F)=0;

proxcost=prox3;
proxcost=proxsigisum;

% proxcost
% 1./prox
% proxcost
% proxcost=proxcost.^2;
% maxproxcost=1e5;
% proxcost(proxcost>maxproxcost)=maxproxcost;
% proxcost
% pause


% prox=1;
% proxt=0;

% visualize for debugging
% for m1=1:F
%     s1=mh(m1).start; e1=mh(m1).end;
%     s1=find(Xo(:,m1),1,'first');    e1=find(Xo(:,m1),1,'last');
%     for m2=m1+1:F
%         s2=mh(m2).start; e2=mh(m2).end;
%         s2=find(Xo(:,m2),1,'first');        e2=find(Xo(:,m2),1,'last');
%
%         mus=max(s1,s2);        mue=min(e1,e2);
%         if mus<=mue
%             evalPerFr=1;
%             tt=mus:mue;
%             %%%% CAREFUL! only for integer tt!!!
%             splinexy1=[Xo(tt,m1)';Yo(tt,m1)'];
%             splinexy2=[Xo(tt,m2)';Yo(tt,m2)'];
%
%             clf; hold on
%             plot3(splinexy1(1,:),splinexy1(2,:),tt,'.-b');
%             plot3(splinexy2(1,:),splinexy2(2,:),tt,'.-r');
%             xlim([ -14070    04981]);ylim([  -14274    01734]); zlim([1 T]);
%             box on; view(3);
%             [m1 m2]
%             tt
%
%
%             d=splinexy1-splinexy2;
%             squares=(d.^2);
%
%             sqrts=sqrt(sum(squares));
%             distalongtraj=sqrts
%             pfl=reshape(proxfromloop(m1,m2,:),1,T)
%             pF*1./prox(m1,m2)
%             prox3(m1,m2)
%             proxsigimax(m1,m2)
%             proxsigisum(m1,m2)
%             pause
%
%             [mindist atframe]=min(distalongtraj);
%
%             %             [m1 m2]
%             %             distalongtraj
% %             prox(m1,m2)=mindist;
% %             proxt(m1,m2)=mus+atframe/evalPerFr-1;
%         end
%     end
% end

% easy vis

if nargout==4
    nMaxAddMerged=50;
    addedMerged=0;
    global used
%     used
    nUsed=length(used);
    merpairs=[];
    
    for m1=1:F
%     for m1=randperm(nUsed)
         mod1=m1;mh1=mh(mod1);

        s1=mh1.start; e1=mh1.end;
        for m2=m1+1:F
%         for  m2=randperm(nUsed)
            mod2=m2; mh2=mh(mod2);
%             [m1 m2]
            if m1==m2, continue; end
            
            if addedMerged>=nMaxAddMerged, break; end            
            
            s2=mh2.start; e2=mh2.end;
            mus=max(s1,s2);        mue=min(e1,e2);
            tt=mus:mue;
            
            [mindist, mintime]=min(proxfromloop(mod1,mod2,:));
            addbuffer=3;        xi=3;
            
%             if m1==35 && m2==53
%             [mod1 mod2]
%             [mintime mus mue mindist]
%         if length(tt)>15
%                         figure(6)
%                         clf
%                         plot(reshape(proxfromloop(m1,m2,:),1,T));
%                         title(sprintf('%d %d %f',m1,m2,mean(proxfromloop(m1,m2,tt))));
%            
% %                         
% %             end
%             mindist
%             mintime
%             [mintime>mus+addbuffer+xi mintime< mue-addbuffer-xi]            
%             mean(proxfromloop(m1,m2,tt)) 
%                          pause
%         end
            if mintime>mus+addbuffer+xi && mintime< mue-addbuffer-xi && mindist<opt.tau*2 && mean(proxfromloop(m1,m2,tt)) > opt.tau 
%                 [m1 m2]
%                 used
%                 [mod1 mod2]
                if    ~isempty(intersect(used,mod1)) && ~isempty(intersect(used,mod2))
                
                    merpairs(end+1,:)=[m1 m2 mus mue mintime];
                    addedMerged=addedMerged+1;
%                         figure(6)
%                         clf
%                         plot(reshape(proxfromloop(m1,m2,:),1,T));
%                         title(sprintf('%d %d %f',m1,m2,mean(proxfromloop(m1,m2,tt))));
%                         pause
                end
                
            end
        end
    end
end
end

function trajTime=findTimeFromSupp(supp, alldpoints,T)
if isempty(supp)
    trajTime=1:T;
else
    %     length(alldpoints.tp)
    %     supp
    supportPtsTs=sort(alldpoints.tp(supp));
    trajTime=[supportPtsTs(1):supportPtsTs(end)];
end
end