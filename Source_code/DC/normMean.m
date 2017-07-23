function [feat featmean stddev]=normMean(feat)
% for track linking

validEntries=find(isfinite(feat));

featmin=min(feat(validEntries));
featmax=max(feat(validEntries));
featmean=mean(feat(validEntries));
featmedian=median(feat(validEntries));
stddev=std(feat(validEntries));
featrange=featmax-featmin;

assert(stddev~=0);


% feat=(feat-featmin)/featrange;
feat=(feat-featmean)/stddev;

end