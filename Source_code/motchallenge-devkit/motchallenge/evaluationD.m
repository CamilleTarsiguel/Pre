%%%% Evaluation Script For Detection %%%%
function [metrics, metricsInfo, additionalInfo] = evaluationD(gtFile, detFile)



% INPUT
% gtFile is the text file containing the ground truth for the sequence
% detFile is the text file containing the results of the detection

gtFile = '/Users/CamilleT/Documents/MATLAB/Pre/groundTruth/SMSQ17short.txt';
detFile = '/Users/CamilleT/Documents/MATLAB/Pre/Source_code/ACF/results/SMSQ17.txt';

% OUTPUT
% metrics contains the following
% [1]   recall	- recall = percentage of detected targets
% [2]   precision	- precision = percentage of correctly detected targets
% [3]   FAR		- number of false alarms per frame
% [4]   GT      - number of GT boxes
% [5]   truepositivies- number of true positives (TP)
% [6]   falsepositives- number of false positives (FP)
% [7]   missed        - number of missed targets (FN)
% [8]	MODA          - N-MODA
% [9]	MODP          - N-MODP


%% Transform files into structures

gtInfo = csvread(gtFile);
stateInfo = csvread(detFile)  ;

%% Evaluation

[metrics, metricsInfo, additionalInfo] = CLEAR_MOD_HUN(gtInfo, stateInfo);