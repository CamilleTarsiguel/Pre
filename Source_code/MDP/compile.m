% --------------------------------------------------------
% MDP Tracking
% Copyright (c) 2015 CVGL Stanford
% Licensed under The MIT License [see LICENSE for details]
% Written by Yu Xiang
% --------------------------------------------------------
%
% compile cpp files
% change the include and lib path if necessary
function compile
OCVRoot = '~/Documents/MATLAB/Pre/Source_code/MDP_Tracking-master/3rd_party/opencv/build';
include = [ '-I', fullfile(OCVRoot, 'include')];
include1 = [ '-I', fullfile(OCVRoot, 'include', 'opencv')];
include2 = [ '-I', fullfile(OCVRoot, 'include','opencv', 'opencv2')];
%lib = ' -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_video';
%eval(['mex lk.cpp -O' include lib]);
mex('lk.cpp', include, include1, include2);
mex distance.cpp 
mex imResampleMex.cpp 
mex warp.cpp

disp('Compilation finished.');
