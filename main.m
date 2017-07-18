%%%% MAIN MAIN MAIN %%%%%
root = cd;
%% Detection only or Tracking?

det = false;
trk = true;

%% Choose dataset

var_MOT16 = false;

var_atrium = false;

var_Smsq1 = false;

var_Smsq2 = false;

var_SMSQ13 = false;

var_SMSQ20 = true;




%% Choose detector

var_peopleDetectorACF = true;

var_CascadeObjectDetector = false;

var_PeopleDetector = false;

var_DPMDetector = false;

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
if var_DPMDetector
    detector = 'DPM';
end


%% Choose tracker

var_MDP = true; % MDP tracking

var_DC = false; % Discrete continuous tracking

var_ROT = false; % Robust online tracking

var_MOT = false; % Baseline with blob analysis

var_MOTv2 = false; % MOT with detectors





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
%% Display mode

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
            
            fr(size1end:size1end + nbIds - 1, 1) =  i * ones( nbIds, 1);
            fr(size1end:size1end + nbIds - 1, 2) = lab';
            fr(size1end:size1end + nbIds - 1, 3) = cpos(1,:)';
            fr(size1end:size1end + nbIds - 1, 4) = cpos(2,:)';
            fr(size1end:size1end + nbIds - 1, 5) = size(1,:)';
            fr(size1end:size1end + nbIds - 1, 6) = size(2,:)';
            fr(size1end:size1end + nbIds - 1, 7) = conf';
            
            size1end = size1end + nbIds;
            
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
end
time = toc(starting);
fprintf("Tracking done in %f \n", toc);
cd(root);