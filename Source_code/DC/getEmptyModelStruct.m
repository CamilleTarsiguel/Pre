function ret=getEmptyModelStruct()
% Trajectory model struct
% similar to a pmpp object with some extensions:
% - start, end:     temporal start and end points
% - labelCost:      corresponding label cost value
% - lcComponents:   label cost components breakdown
% - lastused:       a counter that states how many iterations ago 
%                   this trajectory was used last time
% - index:          piece indexes evaluated at each frame start:end
% - ptindex:        point (detection) indexes that 'belong' to this spline

ret = struct('form',{},'breaks',{},'coefs',{},'pieces',{},'order',{},'dim',{},'bspline',{}, ...
    'start',{},'end',{},'labelCost',{},'lcComponents',{},'lastused',{},'index',{},'ptindex',{},'ptvindex',{});

end