function Dmat = getUnarySpline(nLabels,nPoints,mh,allpoints,opt)
% compute unary energy term for all data points and all labels
% Eq. (7), PAMI


% global opt
dphi=opt.outlierCost;

if isempty(mh)
    Dmat = dphi * ones(nLabels,nPoints);
    return
end

T=opt.seqLength;
normfac=1;
if opt.track3d, normfac=1000; end


Dmat=zeros(nLabels,nPoints);

exthead=2;
exttail=2;
allpind=1:nPoints;

for l=1:nLabels-1
    %
    splstart=mh(l).start-exttail;
    splend=mh(l).end+exthead;
    splineTimespan=max(1,splstart):min(T,splend);
    %
    %         spltime=splstart:splend;
    inSplineTimespan=find(allpoints.tp>=splstart & allpoints.tp<=splend);
%     allothers=setdiff(1:nPoints,inSplineTimespan);
    notinSplineTimespan=true(1,nPoints);
    notinSplineTimespan(inSplineTimespan)=0;
    allothers= allpind & notinSplineTimespan;
    
    %
    pt=allpoints.tp(inSplineTimespan);
    xt=ppval(mh(l),pt)'; %pts on spline
    
    px=allpoints.xp(inSplineTimespan);
    py=allpoints.yp(inSplineTimespan);
    sp=allpoints.sp(inSplineTimespan);
    
    datapts=[px; py]';
    alldists=xt-datapts;
    
    alldists=alldists';
    
    alldists=sqrt(sum(alldists.^2))/normfac; % L2 norm in meter or pixels
    
    
    alln=alldists.^2;

    
    allnlor=-log(1./(alldists.^2+1));

    switch(opt.dataFunction)
        case 1
            Dmat(l,inSplineTimespan)=(alldists .* sp); % L2 norm
        case 2
            Dmat(l,inSplineTimespan)=(alln .* sp); % L2 norm squared
        case 3
            Dmat(l,inSplineTimespan)=(allnlor .* sp);   % Lorentzian
        case 4
            k=opt.conOpt.enParEdat(3);
            alldistsPH = k * (sqrt(1+(alldists./k).^2)-1); % Pseudo-Huber
            alldistsC = sqrt(k+alldists.^2); % Charb. simple
            Dmat(l,inSplineTimespan)=(alldistsC .* sp);   
    end
        
    Dmat(l,allothers)=1e5;
end

uF=opt.unaryFactor;

Dmat=uF*Dmat;


Dmat(nLabels,1:nPoints)=dphi;
Dmat(nLabels,1:nPoints)=dphi*allpoints.sp;
%     Dmat=Dmat+1;
Dmat(Dmat>1e6)=1e6;
Dmat(Dmat<-1e5)=-1e5;

% Dmat=round(Dmat);

end