% --------------------------------------------------------
% MDP Tracking
% Copyright (c) 2015 CVGL Stanford
% Licensed under The MIT License [see LICENSE for details]
% Written by Yu Xiang
% --------------------------------------------------------
%
% build the dres structure for images
function dres_image = read_dres_image_1(opt)

% dres_image.x = zeros(100, 1);
% dres_image.y = zeros(100, 1);
% dres_image.w = zeros(100, 1);
% dres_image.h = zeros(100, 1);
% dres_image.I = cell(100, 1);
% dres_image.Igray = cell(100, 1);
reader = VideoReader(opt.video);

i=1;
while hasFrame(reader)
    
    I = readFrame(reader);

    dres_image.x(i) = 1;
    dres_image.y(i) = 1;
    dres_image.w(i) = size(I, 2);
    dres_image.h(i) = size(I, 1);
    if i ==1
        dres_image.I{i} = I;
        dres_image.Igray{i} = rgb2gray(I);
    end
    
    i = i + 1 ;
end