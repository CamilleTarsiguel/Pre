% --------------------------------------------------------
% MDP Tracking
% Copyright (c) 2015 CVGL Stanford
% Licensed under The MIT License [see LICENSE for details]
% Written by Yu Xiang
% --------------------------------------------------------
%
% test on MOT benchmark
function [dres_track, dres_image, dres_det]  = MOT_test_1(Filename, detector, display)


opt = globals_1(Filename);

global NomDeLaVideo
NomDeLaVideo = Filename;


mot2d_train_seqs = {'TUD-Stadtmitte', 'TUD-Campus', 'PETS09-S2L1', 'ETH-Bahnhof', 'ETH-Sunnyday', 'ETH-Pedcross2'};


% training and testing pairs
%seq_idx_train = {{1, 2}};
seq_idx_train = {1, 2, 3, 4, 5, 6};

N = numel(seq_idx_train);

test_time = 0;

        %idx_train = seq_idx_train{i};
        % load tracker from file
        filename = sprintf('%s/%s_tracker.mat', opt.results, mot2d_train_seqs{1});
        object = load(filename);
        tracker = object.tracker;
        %fprintf('load tracker from file %s\n', filename);
       
    
    % testing
    
    
   
        %fprintf('Testing ...');
        tic;
        [dres_track, dres_image, dres_det] = MDP_test_1(tracker, Filename, detector, display);
        test_time = toc;

%fprintf('Total time for testing: %f\n', test_time);