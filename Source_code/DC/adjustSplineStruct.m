function outspline = ...
    adjustSplineStruct(sfit, sstart, send, allpoints, T, lastused, index, ptindex, ptvindex)
	% we need to append some fields to our spline struct

outspline=sfit;
outspline.start=sstart; outspline.end=send;
[outspline.labelCost, outspline.lcComponents]= ...
    getSplineGoodness(outspline,1,allpoints,T);
outspline.lastused=lastused;
outspline.index=index;
outspline.ptindex=ptindex;
outspline.ptvindex=ptvindex;

end