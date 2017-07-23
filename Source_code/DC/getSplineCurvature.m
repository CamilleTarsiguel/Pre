function curvature = getSplineCurvature(spl, tt, opt)
% Compute the total curvature of a spline sampled at
% discrete time steps
%
% if explicit timesteps are not given,
% the entire spline is considered
%
% see Wikiepedia for definition
% http://en.wikipedia.org/wiki/Curvature#Curvature_of_space_curves
% 

if nargin<2
    tt=spl.start:spl.end;
end


s1=ppdiff(spl,1);s2=ppdiff(spl,2);

samplrate=1;
% tt=linspace(tt(1),tt(end),1/samplrate*length(tt));

A1 = ppval(s1,tt); % first derivative
A2 = ppval(s2,tt); % second derivative
A1x=A1(1,:);A1y=A1(2,:);
A2x=A2(1,:);A2y=A2(2,:);
% curvature=sqrt((-A2y).^2 + A2x.^2 + (A2y.*A1x - A2x.*A1y).^2)./((A1x.^2 + A1y.^2 + 1).^(3/2));
% curvature=abs(A1x.*A2y - A1y.*A2x)./((A1x.^2+A1y.^2).^(3/2));
%  curvature=abs(A1x.*A2y - A1y.*A2x)./(A1x.^2+A1y.^2) * samplrate;


% pseudo Huber
%  k=.1;
%  nomval=A1x.*A2y - A1y.*A2x;
%  nomval=k^2 * (sqrt(1+(nomval./k).^2)-1);
% curvature= nomval ./ (A1x.^2+A1y.^2) * samplrate;
% sum(curvature)

% Simple Charbonnier
k=1; epsil=.1;
nomval=A1x.*A2y - A1y.*A2x;
nomval= sqrt(epsil+(nomval).^2);
% curvature= nomval ./ (A1x.^2+A1y.^2) * samplrate;
curvature= nomval ./ (A1x.^2+A1y.^2 + epsil) * samplrate;
end