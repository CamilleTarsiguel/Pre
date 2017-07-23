function checkSplines(spl)
% Check if all splines are cubic

if any([spl(:).order]~=4), error('All splines must be cubic'); end