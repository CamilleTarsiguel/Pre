%%%% Script débug MDP %%%%%%


%%%% MAIN MAIN MAIN %%%%%


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


%% Loading video

if var_MOT16 % video in Shanghai
    filename = '~/Documents/MATLAB/Pre/data/MOT16-04.mp4';
end
if var_Smsq1 % Video far away
    filename = '~/Documents/MATLAB/Pre/data/SmSq1.m4v';
end
if var_SMSQ13 % Video from the Kitta
    filename = '~/Documents/MATLAB/Pre/data/SMSQ13.MOV';
end
if var_SMSQ20 % Video from the Kitta (group of people)
    filename = '~/Documents/MATLAB/Pre/data/SMSQ20short.avi';
end

%% Display mode

display = false;
starting = tic;
%% Tracking ( + detection for online)
fprintf('Tracking Start\n');
    fprintf('MDP Tracking ... \n');
    [dres_track, dres_image, dres_det] = MOT_test_1(filename, detector);
time = toc(starting);
fprintf("Tracking done in %f \n", toc);
