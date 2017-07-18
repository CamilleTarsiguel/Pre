%% Train Faster R-CNN Vehicle Detector

%%
% Load training data.
data = load('fasterRCNNVehicleTrainingData.mat');
trainingData = load('/Users/CamilleT/Documents/MATLAB/Pre/Source_code/7/Det/det.mat');

trainingData.imageFilename = fullfile(toolboxdir('vision'),'visiondata', ...
    trainingData.imageFilename);

%%
% Setup network layers.
layers = data.layers

%%
% Configure training options.
%
% * Lower the InitialLearningRate to reduce the rate at which network
%   parameters are changed.
% * Set the CheckpointPath to save detector checkpoints to a temporary
%   directory. Change this to another location if required.
% * Set MaxEpochs to 1 to reduce example training time. Increase this
%   to 10 for proper training.
options = trainingOptions('sgdm', ...
    'InitialLearnRate', 1e-6, ...
    'MaxEpochs', 1, ...
    'CheckpointPath', tempdir);

%%
% Train detector. Training will take a few minutes.
detector = trainFasterRCNNObjectDetector(trainingData, layers, options)

%%
% Test the Fast R-CNN detector on a test image.
img = imread('highway.png');

%%
% Run detector.
[bbox, score, label] = detect(detector, img);

%%
% Display detection results.
detectedImg = insertShape(img, 'Rectangle', bbox);
figure
imshow(detectedImg)