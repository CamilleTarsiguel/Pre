function [prox, proxt, proxcost]= ...
    getOneSplineProximity(mh,T,stInfo,m,opt)
% compute distances from one spline (m) to all others
%


F=length(mh);


% global opt sceneInfo
pF=opt.proxcostFactor;
proxcost=zeros(1,F);
prox=Inf*ones(1,F);
proxt=zeros(1,F);
if ~pF
    return;
end


% pairwise distance between splines

tt=1:T;

% compute all spline positions
if nargin<3,    stInfo=[]; end

if isempty(stInfo)
    X=zeros(T,F);
    Y=zeros(T,F);
    for m1=1:F
        tt=mh(m1).start:mh(m1).end;
        splinexy=ppval(mh(m1),tt);
        X(tt,m1)=splinexy(1,:)';
        Y(tt,m1)=splinexy(2,:)';
    end
else
    X=stInfo.X; Y=stInfo.Y;
end

% vectorized version
X(~X)=Inf; Y(~Y)=Inf;

% prox=Inf*ones(F,T);

tmp=repmat(X(:,m),1,F);
dx=(tmp-X).^2;
tmp=repmat(Y(:,m),1,F);
dy=(tmp-Y).^2;
prox=sqrt(dx+dy);
    % inf-inf = NaN
    prox(isnan(prox))=Inf;

proxfirst=prox;
[prox proxt]=min(prox);

prox3=pF*sum(1./proxfirst);
% for t=1:T
%     tmp=repmat(X(t,:),F,1);
%     ax=tmp(:)';
%     tmp=tmp';
%     bx=tmp(:)';
%     
%     tmp=repmat(Y(t,:),F,1);
%     ay=tmp(:)';
%     tmp=tmp';
%     by=tmp(:)';
%     
%     dx=(ax-bx).^2;
%     dy=(ay-by).^2;
%     d=sqrt(dx+dy);
%     
%     prox2(:,:,t)=reshape(d,F,F);    
% end
% [prox2 proxt2]=min(prox2,[],3);
% prox2=triu(prox2); prox2(isinf(prox2))=Inf;
% prox2(~~(tril(ones(F))))=Inf;
% prox
% prox2
% prox=prox2;
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


proxcost=pF*1./prox;
proxcost=prox3;

siga=opt.conOpt.enParEexc(1);
sigb=(opt.conOpt.enParEexc(2)/2)*siga;
sigscale=pF;
proxsigi=sigscale-sigscale*1./(1+exp(-siga*proxfirst+sigb));
proxsigimax=max(proxsigi);
proxsigisum=sum(proxsigi);

proxsigimax(m)=0;
proxsigisum(m)=0;
% m
% proxfirst
% proxsigisum
% proxsigimax
% prox3
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