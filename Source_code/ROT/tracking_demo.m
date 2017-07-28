%% Robust Online Multi-Object Tracking based on Tracklet Confidence
%% and Online Discriminative Appearance Learning (CVPR2014)
% Last updated date: 2014. 07. 27
%% Copyright (C) 2014 Seung-Hwan Bae
%% All rights reserved.

function [Trk_sets, detections, all_mot] = tracking_demo(filename, detector, display)
%clc;
%clear all
base = [pwd, '/'];
addpath(genpath(base));
fprintf('You choose the video %s \n', filename);
global Filename;
Filename = filename;
mot_setting_params;

% Setting the detector
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


 % Choose parameters

maxH = 100;
minH = 2;
maxW = 70;
minW = 2;
 
% Loading Video

reader = VideoReader(video_path);
frame_end = reader.NumberOfFrames;
reader = VideoReader(video_path);


fprintf('Finished changing struture \n');

% 1:ILDA, 0: No-ILDA (faster)
% To use ILDA, refer to README.
param.use_ILDA = 0;


tic;
frame_start = 1;

All_Eval = [];
cct = 0;
Trk = []; Trk_sets = []; all_mot =[];

%% Initialization Tracklet
fprintf('Initialization Tracklets \n');
tstart1 = tic;
init_frame = frame_start + param.show_scan;

for i=1:init_frame
    frame = readFrame(reader);
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
    detections(i).x = bboxes(:,1);
    detections(i).y = bboxes(:,2);
    detections(i).w = bboxes(:,3);
    detections(i).h = bboxes(:,4);
    Obs_grap(i).iso_idx = ...
        ones(size(detections(i).x));
    Obs_grap(i).child = [];
    Obs_grap(i).iso_child =[];
    init_img_set{i} = frame;

end

[Obs_grap] = mot_pre_association(detections,Obs_grap,frame_start,init_frame);



[Trk,param,Obs_grap] = MOT_Initialization_Tracklets(init_img_set,Trk,detections,param,...
    Obs_grap,init_frame);
fprintf('Initialization Tracklets finished \n');
%% Tracking
disp('Tracking objects...');
for fr = init_frame+1:frame_end
    if fr ==10
        param.new_thr = 5;
    end
    rgbimg = readFrame(reader);
    
    % Detection
    if var_peopleDetectorACF
        [bboxes, score ] = detect(Detector, rgbimg);
    end
    if var_CascadeObjectDetector
        bboxes = step(Detector, rgbimg);
        score = ones(size(bboxes,1),1);
    end
    
    if var_PeopleDetector
        [bboxes, score ] = step(Detector, rgbimg);
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
    detections(fr).x = bboxes(:,1);
    detections(fr).y = bboxes(:,2);
    detections(fr).w = bboxes(:,3);
    detections(fr).h = bboxes(:,4);
    init_img_set{fr} = rgbimg;
    
    %% Local Association
    [Trk, Obs_grap, Obs_info] = MOT_Local_Association(Trk, detections, Obs_grap, param, ILDA, fr, rgbimg);
    %% Global Association
    [Trk, Obs_grap] = MOT_Global_Association(Trk, Obs_grap, Obs_info, param, ILDA, fr);
    %% Tracklet Confidence Update
    [Trk] = MOT_Confidence_Update(Trk,param,fr, param.lambda);
    [Trk] = MOT_Type_Update(rgbimg,Trk,param.type_thr,fr);
    %% Tracklet State Update & Tracklet Model Update
    [Trk] = MOT_State_Update(Trk, param, fr);
    %% New Tracklet Generation
    [Trk, param, Obs_grap] = MOT_Generation_Tracklets(init_img_set,Trk,detections,param,...
        Obs_grap,fr);
    
    %% Incremental subspace learning
    %     if param.use_ILDA
    %         [ILDA] = MOT_Online_Appearance_Learning(rgbimg, img_path, img_List, fr, Trk, param, ILDA);
    %     end
    
    %% Tracking Results
    [Trk_sets] = MOT_Tracking_Results(Trk,Trk_sets,fr);
    disp([sprintf('Tracking:Frame_%04d',fr)]);
end
%%
toc
disp('Tracking done...');
TotalTime = toc(tstart1);
AverageTime = TotalTime/(frame_start + frame_end);

%% Draw Tracking Results
out_path = 'Results/SMSQ17/';

DrawOption.isdraw = display;
DrawOption.iswrite = 1;
DrawOption.new_thr = param.new_thr;

% Box colors indicate the confidences of tracked objects
% High (Red)-> Low (Blue)
[all_mot] = MOT_Draw_Tracking(Trk_sets, out_path, video_path, DrawOption);
close all;
disp([sprintf('Average running time:%.3f(sec/frame)', AverageTime)]);


end



