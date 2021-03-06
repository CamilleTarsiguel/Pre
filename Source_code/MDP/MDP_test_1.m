% --------------------------------------------------------
% MDP Tracking
% Copyright (c) 2015 CVGL Stanford
% Licensed under The MIT License [see LICENSE for details]
% Written by Yu Xiang
% --------------------------------------------------------
%
% testing MDP
function [dres_track, dres_image, dres_det] = MDP_test_1(tracker,Filename, detector, display)

is_show = display;   % set is_show to 1 to show tracking results in testing
is_save = 0;   % set is_save to 1 to save tracking result
is_text = 0;   % set is_text to 1 to display detailed info
is_pause = 0;  % set is_pause to 1 to debug

opt = globals_1(Filename);
opt.is_text = is_text;
opt.exit_threshold = 0.7;

if is_show
    close all;
end

seq_name = extractBefore(extractAfter(Filename, 'data/'),'.');

%% build the dres structure for images
filename = sprintf('%s/%s_dres_image.mat', opt.results, seq_name);
dres_image = read_dres_image_1(opt);
%fprintf('read images done\n');
seq_num = length(dres_image.x);
% % read detections
% filename = opt.det;
% dres_det = read_mot2dres(filename);
%
reader = VideoReader(opt.video);
%readFrame(reader);
%% Choosing the detector and initializing the detector
switch detector
    case 'ACF'
    var_peopleDetectorACF = true;
    var_CascadeObjectDetector = false;
    var_PeopleDetector = false;
    case 'Cascade'
    var_peopleDetectorACF = false;
    var_CascadeObjectDetector = true;
    var_PeopleDetector = false;
    case 'HOG'
    var_peopleDetectorACF = false;
    var_CascadeObjectDetector = false;
    var_PeopleDetector = true;
    otherwise
        error('This detector is not supported');
end
if var_peopleDetectorACF
    Detector = peopleDetectorACF('caltech-50x21');
end
if var_CascadeObjectDetector
    Detector = vision.CascadeObjectDetector('UpperBody');
end
if var_PeopleDetector
    Detector = vision.PeopleDetector('UprightPeople_96x48');
end


%% load the trained model
if nargin < 3
    object = load('results/TUD-Campus_tracker.mat');
    tracker = object.tracker;
end

%% intialize tracker
I = dres_image.I{1};
tracker = MDP_initialize_test(tracker, size(I,2), size(I,1), is_show);

%% for each frame
trackers = [];
id = 0;
for fr = 1:seq_num
    if fr ~= 1
        dres_image.I{fr-1} = [];
        dres_image.Igray{fr-1} = [];
    end
    dres_image.I{fr} = readFrame(reader);
    dres_image.Igray{fr} = rgb2gray(dres_image.I{fr});
    frame = dres_image.I{fr};
    if is_text
        fprintf('frame %d\n', fr);
    else
        fprintf('.');
        if mod(fr, 100) == 0
            fprintf('\n');
        end        
    end
    [score, bboxes] = detectObjects(Detector, frame, var_peopleDetectorACF, var_CascadeObjectDetector, var_PeopleDetector);
    dres_det.fr = fr * ones(size(bboxes, 1),1);
    dres_det.id = -1 * ones(size(bboxes, 1),1);
    dres_det.x = bboxes(:,1);
    dres_det.y = bboxes(:,2);
    dres_det.w = bboxes(:,3);
    dres_det.h = bboxes(:,4);
    dres_det.r = score;
    
    dres = dres_det;
    
    % nms
%     boxes = [dres.x dres.y dres.x+dres.w dres.y+dres.h dres.r];
%     index = nms_new(boxes, 0.6);
%     dres = sub(dres, index);
    
    dres = MDP_crop_image_box(dres, dres_image.Igray{fr}, tracker);
    
    if is_show
        figure(1);
        
        

        % show detections
        subplot(2, 2, 2);
        show_dres(fr, dres_image.I{fr}, 'Detections', dres_det);
    end
%     
    % sort trackers
    index_track = sort_trackers(trackers);
    
    % process trackers
    for i = 1:numel(index_track)
        ind = index_track(i);
        
        if trackers{ind}.state == 2
            % track target
            trackers{ind} = track(fr, dres_image, dres, trackers{ind}, opt);
            % connect target
            if trackers{ind}.state == 3
                [dres_tmp, index] = generate_initial_index(trackers(index_track(1:i-1)), dres);
                dres_associate = sub(dres_tmp, index);
                trackers{ind} = associate(fr, dres_image,  dres_associate, trackers{ind}, opt);
            end
        elseif trackers{ind}.state == 3
            % associate target
            [dres_tmp, index] = generate_initial_index(trackers(index_track(1:i-1)), dres);
            dres_associate = sub(dres_tmp, index);    
            trackers{ind} = associate(fr, dres_image, dres_associate, trackers{ind}, opt);
        end
    end
    
    % find detections for initialization
    [dres, index] = generate_initial_index(trackers, dres);
    for i = 1:numel(index)
        % extract features
        dres_one = sub(dres, index(i));
        f = MDP_feature_active(tracker, dres_one);
        % prediction
        label = svmpredict(1, f, tracker.w_active, '-q');
        % make a decision
        if label < 0
            continue;
        end
        
        % reset tracker
        tracker.prev_state = 1;
        tracker.state = 1;            
        id = id + 1;
        
        trackers{end+1} = initialize(fr, dres_image, id, dres, index(i), tracker);
    end
    
    % resolve tracker conflict
    trackers = resolve(trackers, dres, opt);    
    
    dres_track = generate_results(trackers);
    if is_show
        figure(1);

        %show tracking results
        subplot(2, 2, 3);
        show_dres(fr, dres_image.I{fr}, 'Tracking', dres_track, 2);

        % show lost targets
        subplot(2, 2, 4);
        show_dres(fr, dres_image.I{fr}, 'Lost', dres_track, 3);

        if is_pause
            pause();
        else
            pause(0.01);
        end
    end  
end

%% write tracking results
filename = sprintf('%s/%s.txt', opt.results, seq_name);
%fprintf('write results: %s\n', filename);
write_tracking_results(filename, dres_track, opt.tracked);
%% Detection

  function [score, bboxes] = detectObjects(Detector, frame,var_peopleDetectorACF, var_CascadeObjectDetector, var_PeopleDetector )

      
      % some parameters
maxH = 100;
minH = 2;
maxW = 70;
minW = 2;

      
        if var_peopleDetectorACF
            [bboxes, score ] = detect(Detector, frame);
        end
        if var_CascadeObjectDetector
            bboxes = step(Detector, frame);
            score = ones(size(bboxes,1),1);
        end
        
        if var_PeopleDetector
            [bboxes, score ] = step(Detector, frame);
        end
        index = [bboxes(:,3)<maxW];
        bboxes = bboxes(index, :);
        score = score(index, 1);
        index = [bboxes(:,3)>minW];
        bboxes = bboxes(index, :);
        score = score(index, 1);
        index = [bboxes(:,4)<maxH];
        bboxes = bboxes(index, :);
        score = score(index, 1);
        index = [bboxes(:,4)>minH];
        bboxes = bboxes(index, :);
        score = score(index, 1);


 

%% sort trackers according to number of tracked frames
function index = sort_trackers(trackers)

sep = 10;
num = numel(trackers);
len = zeros(num, 1);
state = zeros(num, 1);
for i = 1:num
    len(i) = trackers{i}.streak_tracked;
    state(i) = trackers{i}.state;
end

index1 = find(len > sep);
[~, ind] = sort(state(index1));
index1 = index1(ind);

index2 = find(len <= sep);
[~, ind] = sort(state(index2));
index2 = index2(ind);
index = [index1; index2];


%% initialize a tracker
% dres: detections
function tracker = initialize(fr, dres_image, id, dres, ind, tracker)

if tracker.state ~= 1
    return;
else  % active

    % initialize the LK tracker
    tracker = LK_initialize(tracker, fr, id, dres, ind, dres_image);
    tracker.state = 2;
    tracker.streak_occluded = 0;
    tracker.streak_tracked = 0;

    % build the dres structure
    dres_one.fr = dres.fr(ind);
    dres_one.id = tracker.target_id;
    dres_one.x = dres.x(ind);
    dres_one.y = dres.y(ind);
    dres_one.w = dres.w(ind);
    dres_one.h = dres.h(ind);
    dres_one.r = dres.r(ind);
    dres_one.state = tracker.state;
    tracker.dres = dres_one;
end


%% track a target
function tracker = track(fr, dres_image, dres, tracker, opt)

% tracked    
if tracker.state == 2
    tracker.streak_occluded = 0;
    tracker.streak_tracked = tracker.streak_tracked + 1;
    tracker = MDP_value(tracker, fr, dres_image, dres, []);

    % check if target outside image
    [~, ov] = calc_overlap(tracker.dres, numel(tracker.dres.fr), dres_image, fr);
    if ov < opt.exit_threshold
        if opt.is_text
            fprintf('target outside image by checking boarders\n');
        end
        tracker.state = 0;
    end    
end



%% associate a lost target
function tracker = associate(fr, dres_image, dres_associate, tracker, opt)

% occluded
if tracker.state == 3
    tracker.streak_occluded = tracker.streak_occluded + 1;
    % find a set of detections for association
    [dres_associate, index_det] = generate_association_index(tracker, fr, dres_associate);
    tracker = MDP_value(tracker, fr, dres_image, dres_associate, index_det);
    if tracker.state == 2
        tracker.streak_occluded = 0;
    end

    if tracker.streak_occluded > opt.max_occlusion
        tracker.state = 0;
        if opt.is_text
            fprintf('target %d exits due to long time occlusion\n', tracker.target_id);
        end
    end
    
    % check if target outside image
    [~, ov] = calc_overlap(tracker.dres, numel(tracker.dres.fr), dres_image, fr);
    if ov < opt.exit_threshold
        if opt.is_text
            fprintf('target outside image by checking boarders\n');
        end
        tracker.state = 0;
    end    
end


%% resolve conflict between trackers
function trackers = resolve(trackers, dres_det, opt)

% collect dres from trackers
dres_track = [];
for i = 1:numel(trackers)
    tracker = trackers{i};
    dres = sub(tracker.dres, numel(tracker.dres.fr));
    
    if tracker.state == 2
        if isempty(dres_track)
            dres_track = dres;
        else
            dres_track = concatenate_dres(dres_track, dres);
        end
    end
end   

% compute overlaps
num_det = numel(dres_det.fr);
if isempty(dres_track)
    num_track = 0;
else
    num_track = numel(dres_track.fr);
end

flag = zeros(num_track, 1);
for i = 1:num_track
    [~, o] = calc_overlap(dres_track, i, dres_track, 1:num_track);
    o(i) = 0;
    o(flag == 1) = 0;
    [mo, ind] = max(o);
    if mo > opt.overlap_sup        
        num1 = trackers{dres_track.id(i)}.streak_tracked;
        num2 = trackers{dres_track.id(ind)}.streak_tracked;
        o1 = max(calc_overlap(dres_track, i, dres_det, 1:num_det));
        o2 = max(calc_overlap(dres_track, ind, dres_det, 1:num_det));
        
        if num1 > num2
            sup = ind;
        elseif num1 < num2
            sup = i;
        else
            if o1 > o2
                sup = ind;
            else
                sup = i;
            end
        end
        
        trackers{dres_track.id(sup)}.state = 3;
        trackers{dres_track.id(sup)}.dres.state(end) = 3;
        if opt.is_text
            fprintf('target %d suppressed\n', dres_track.id(sup));
        end
        flag(sup) = 1;
    end
end


