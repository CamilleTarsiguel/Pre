%%%% Evaluation Script For Tracking %%%%
function [metrics, metricsInfo, additionalInfo] = evaluationT(gtFile, trkFile)

%% Set files

%trkFile = '/Users/CamilleT/Documents/Codes_annexes/MDP_Tracking/results/SMSQ20.txt';
%trkFile = '/Users/CamilleT/Documents/MATLAB/Pre/Source_code/KLT/results/SMSQ17.txt';
%gtFile = '/Users/CamilleT/Documents/MATLAB/Pre/groundTruth/SMSQ20short.txt';
%gtFile = '/Users/CamilleT/Documents/MATLAB/Pre/Source_code/ROT/results/SMSQ20short.txt';

%% Set parameters

options.eval2d = 1;
options.eval3d = 0;
%options.td = .5; % threshold

%% Transform files into structures

gtInfo = convertTXTToStruct(gtFile);
stateInfo = convertTXTToStruct(trkFile);

%% Evaluation

[metrics, metricsInfo, additionalInfo] = CLEAR_MOT_HUN(gtInfo, stateInfo, options);
end