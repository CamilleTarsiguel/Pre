% --------------------------------------------------------
% MDP Tracking
% Copyright (c) 2015 CVGL Stanford
% Licensed under The MIT License [see LICENSE for details]
% Written by Yu Xiang
% --------------------------------------------------------
%
% initialization for testing
function tracker = MDP_initialize_test(tracker, image_width, image_height, is_show)

% normalization factor for features
tracker.image_width = image_width;
tracker.image_height = image_height;
tracker.max_width = 50;
tracker.max_height = 100;
tracker.max_score = 166.4;

tracker.streak_tracked = 0;
tracker.is_show = is_show;