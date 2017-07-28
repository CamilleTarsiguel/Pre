% --------------------------------------------------------
% MDP Tracking
% Copyright (c) 2015 CVGL Stanford
% Licensed under The MIT License [see LICENSE for details]
% Written by Yu Xiang
% --------------------------------------------------------
%
% test on MOT benchmark
function MOT_test_1

opt = globals_1();

mot2d_train_seqs = {'TUD-Stadtmitte', 'TUD-Campus'};


% training and testing pairs
seq_idx_train = {{1, 2}};

N = numel(seq_idx_train);

test_time = 0;
for i = 1:N
        idx_train = seq_idx_train{i};
        % load tracker from file
        filename = sprintf('%s/%s_tracker.mat', opt.results, mot2d_train_seqs{idx_train{end}});
        object = load(filename);
        tracker = object.tracker;
        fprintf('load tracker from file %s\n', filename);
       
    
    % testing
    
    % number of testing sequences
   
        fprintf('Testing ...');
        tic;
        MDP_test_1(tracker);
        test_time = test_time + toc;
end

fprintf('Total time for testing: %f\n', test_time);