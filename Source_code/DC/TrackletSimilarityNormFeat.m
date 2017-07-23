%% ETH only
simmatfile=sprintf('tmp/linked/simMat-clnr%03d-s%d.mat',clusternr,51);
load(simmatfile);

nrmlzDist=sumdist.^2;nrmlzScale=sumdistS.^2;
nrmlzBhat=distBhatt;nrmlzGaps=gaps;


D1=nrmlzDist(:);S1=nrmlzScale(:);
B1=nrmlzBhat(:);G1=nrmlzGaps(:);

simmatfile=sprintf('tmp/linked/simMat-clnr%03d-s%d.mat',clusternr,53);
load(simmatfile);

nrmlzDist=sumdist.^2;nrmlzScale=sumdistS.^2;
nrmlzBhat=distBhatt;nrmlzGaps=gaps;


D2=nrmlzDist(:);S2=nrmlzScale(:);
B2=nrmlzBhat(:);G2=nrmlzGaps(:);


D=[D1;D2];S=[S1;S2];
B=[B1;B2];G=[G1;G2];


[min(D(isfinite(D))) max(D(isfinite(D)))]
[D featmeanD stddevD]=normMean(D);
[S featmeanS stddevS]=normMean(S);
[B featmeanB stddevB]=normMean(B);
[G featmeanG stddevG]=normMean(G);
[min(D(isfinite(D))) max(D(isfinite(D)))]

save(sprintf('tmp/linked/simMat-clnr%03d-featnorm.mat',clusternr),'featmean*','stddev*');
