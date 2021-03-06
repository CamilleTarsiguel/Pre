%%% Evaluation For Detection And / Or tracking %%%

root = cd;


%% Detection or Tracking ?

det = true;
trk = false;




%% Choose dataset

var_MOT16 = false;

var_atrium = false;

var_Smsq1 = false;

var_Smsq2 = false;

var_SMSQ13 = false;

var_SMSQ20 = true;

var_SMSQ24 = false;

var_SMSQ16 = false;

var_SMSQ14 = false;

var_SMSQ17 = false;



% Choose detector

var_peopleDetectorACF = true;

var_CascadeObjectDetector = false;

var_PeopleDetector = false;

% make the variable

if var_peopleDetectorACF
    detector = 'ACF';
end
if var_CascadeObjectDetector
    detector = 'Cascade';
end
if var_PeopleDetector
    detector = 'HOG';
end




%% Choose tracker

var_MDP = false; % MDP tracking

var_DC = false; % Discrete continuous tracking

var_ROT = true; % Robust online tracking

var_MOT = false; % Baseline with blob analysis

var_MOTv2 = false; % MOT with detectors

var_KLT = false; % Online tracker based on KLT





%% Loading video

if var_MOT16 % video in Shanghai
    filename = fullfile(root,'data/MOT16-04.mp4');
    videoname = 'MOT16-04';
end
if var_Smsq1 % Video far away
    filename = fullfile(root,'/data/SmSq1.m4v');
    videoname = 'SmSq1';
end
if var_SMSQ13 % Video from the Kitta
    filename = fullfile(root,'data/SMSQ13.MOV');
    videoname = 'SMSQ13';
end
if var_SMSQ20 % Video from the Kitta (group of people)
    filename = fullfile(root,'data/SMSQ20short50.avi');
    videoname = 'SMSQ20short';
end
if var_SMSQ24 % Video from the Kitta (group of people)
    filename = fullfile(root,'data/SMSQ24.MOV');
    videoname = 'SMSQ24';
end
if var_SMSQ16 % Video from the Kitta (group of people)
    filename = fullfile(root,'data/SMSQ16.MOV');
    videoname = 'SMSQ16';
end
if var_SMSQ14 % Video from the Kitta (group of people)
    filename = fullfile(root,'data/SMSQ14.MOV');
    videoname = 'SMSQ14';
end
if var_SMSQ17 % Video from the Kitta (group of people)
    filename = fullfile(root,'data/SMSQ17short.avi');
    videoname = 'SMSQ17short';
end


%% Do you need to run the code now?

Y = false; % Yes I need to
N = true;  % No it's okay



%% You said yes
if Y
% Display mode

display = false;
starting = tic;
%% Detection (offline)
if (var_DC||det) % For now the detectors are still not turned into online trackers
    detect_people(filename,videoname, detector, display);
end

%% Tracking ( + detection for online)
if trk
    fprintf('Tracking Start\n');
    if var_MDP
        cd('Source_code/MDP/');
        fprintf('MDP Tracking ... \n');
        [dres_track, dres_image, dres_det] = MOT_test_1(filename, detector, display);
    end
    
    if var_DC
        cd('Source_code/DC');
        fprintf('DC Tracking ... \n');
        % Execute
        [metrics2d, metrics3d, allen, stateInfo, sceneInfo]=swDCTracker('scenes/perso.ini','config/default2dSimple.ini');
        
        
        % Save results
        numframes = stateInfo.F;
        fr = [];
        nbIds = size( stateInfo.X, 2);
        size1end = 1;
        
        for i=1:numframes
            fr(size1end:size1end + nbIds - 1, 1) =  i * ones( nbIds, 1);
            fr(size1end:size1end + nbIds - 1, 2) = (1:nbIds)';
            fr(size1end:size1end + nbIds - 1, 3) = stateInfo.Xi(i,:)';
            fr(size1end:size1end + nbIds - 1, 4) = stateInfo.Yi(i,:)';
            fr(size1end:size1end + nbIds - 1, 5) = stateInfo.W(i,:)';
            fr(size1end:size1end + nbIds - 1, 6) = stateInfo.H(i,:)';
            size1end = size1end + nbIds;
        end
        
        csvwrite(strcat('results/',videoname, '.txt'), fr);
        
    end
    if var_ROT
        cd('Source_code/ROT');
        fprintf('ROT Tracking ... \n');
        [Trk_sets, detections, all_mot] = tracking_demo(filename, detector, display);
        
        numframes = length(Trk_sets);
        fr = [];
        size1end = 1;
        for i = 1 : numframes
            
            cpos = cell2mat(all_mot.cpos(i));
            size = cell2mat(all_mot.size(i));
            lab = cell2mat(all_mot.lab(i));
            conf = cell2mat(all_mot.conf(i));
            nbIds = length(cpos);
            if ~sempty(conf)
            fr(size1end:size1end + nbIds - 1, 1) =  i * ones( nbIds, 1);
            fr(size1end:size1end + nbIds - 1, 2) = lab';
            fr(size1end:size1end + nbIds - 1, 3) = cpos(1,:)';
            fr(size1end:size1end + nbIds - 1, 4) = cpos(2,:)';
            fr(size1end:size1end + nbIds - 1, 5) = size(1,:)';
            fr(size1end:size1end + nbIds - 1, 6) = size(2,:)';
            fr(size1end:size1end + nbIds - 1, 7) = conf';
            
            size1end = size1end + nbIds;
            end
        end
        csvwrite(strcat('results/',videoname, '.txt'), fr);
        
    end
    if var_MOT
        memory = MOT(filename, display);
        csvwrite(fullfile('Source_code', 'MOT','results', strcat(videoname, '.txt')),memory);
    end
    if var_MOTv2
        memory2 = MOT_v2(filename, detector, display);
        csvwrite(fullfile('Source_code', 'MOT_v2','results', strcat(videoname, '.txt')),memory2);
    end
    if var_KLT
        memory3 = MOT_KLT(filename, display);
        end1 = 1;
        for i = 1 : size(memory3,2) % id
            l = size(memory3(i).bboxes,1);
            fr(end1:end1+l-1,2)= i * ones(l,1);
            fr(end1:end1+l-1,1)= memory3(i).frame;
            fr(end1:end1+l-1,3:6)= memory3(i).bboxes;
            fr(end1:end1+l-1,7)= 1 * ones(l,1);
            
            end1 = end1 + l;
        end
        csvwrite(fullfile('Source_code', 'KLT','results', strcat(videoname, '.txt')),fr);
    end
end
time = toc(starting);
fprintf("Tracking done in %f \n", time);
reader = VideoReader(filename);
nbf = reader.NumberOfFrames;
speed = nbf/time;
fprintf("Speed is %d fps \n", speed);
cd(root);
end


%% You said no



if det
    if var_peopleDetectorACF
        detFile = fullfile(root,'Source_code', 'ACF', 'results', strcat(videoname ,'.txt'));
    end
    if var_CascadeObjectDetector
        detFile = fullfile(root,'Source_code', 'Cascade', 'results', strcat(videoname ,'.txt'));
    end
    if var_PeopleDetector
        detFile = fullfile(root,'Source_code', 'HOG', 'results', strcat(videoname ,'.txt'));
    end
    gtFile = fullfile(root, 'groundTruth',strcat(videoname, '.txt'));
    
    cd(fullfile(root,'Source_code/motchallenge-devkit/motchallenge/'))
    addpath(fullfile(root,'Source_code/motchallenge-devkit/motchallenge/', 'utils'));
    [metrics, metricsInfo, additionalInfo] = evaluationD(gtFile, detFile)
end
if trk
    
    if var_MDP 
        trkFile = fullfile(root, 'Source_code', 'MDP', 'results',strcat(videoname, '.txt')) ;

    end
    
    if var_DC 
    end
    
    if var_ROT
        trkFile = fullfile(root, 'Source_code', 'ROT', 'results',strcat(videoname, '.txt'));
    end

    if var_MOT
        trkFile = fullfile(root, 'Source_code', 'MOT', 'results',strcat(videoname, '.txt'));
    end

    if var_MOTv2
        trkFile = fullfile(root, 'Source_code', 'MOT_v2', 'results',strcat(videoname, '.txt'));
    end
    if var_KLT
        trkFile = fullfile(root, 'Source_code', 'KLT', 'results',strcat(videoname, '.txt'));
    end
    gtFile = fullfile(root, 'groundTruth',strcat(videoname, '.txt'));
    %gtFile = fullfile(root, 'Source_code', 'ROT', 'results',strcat(videoname, '.txt'));
    cd(fullfile(root,'Source_code/motchallenge-devkit/motchallenge/'))
    addpath(fullfile(root,'Source_code/motchallenge-devkit/motchallenge/','utils'));
    [metrics, metricsInfo, additionalInfo] = evaluationT(gtFile, trkFile);
    metrics
    
end
    
cd(root);