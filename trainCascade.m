%% Train Cascade Object detector
% Load the positive samples data from a MAT file. The file contains
% a table specifying bounding boxes for several object categories.
% The table was exported from the Training Image Labeler app.
%%
% Load positive samples.
%load('stopSignsAndCars.mat');
%%
% Select the bounding boxes for stop signs from the table.
positiveInstances = trainingTable;
%%
% Add the image directory to the MATLAB path.
imDir = '/Users/CamilleT/Documents/MATLAB/Pre/Seq_2_AVI/data/Caltech/train04/images/';
addpath(imDir);
%%
% Specify the foler for negative images.
negativeFolder = '/Users/CamilleT/Documents/MATLAB/Pre/Seq_2_AVI/data/Caltech/train04/neg/';
%%
% Create an |imageDatastore| object containing negative images.
negativeImages = imageDatastore(negativeFolder);
%%
% Train a cascade object detector called 'stopSignDetector.xml'
% using HOG features.
% NOTE: The command can take several minutes to run.
tic;
%trainCascadeObjectDetector('CascadePeople.xml', 'resume');
 trainCascadeObjectDetector('CascadePeople.xml',positiveInstances, ...
     negativeFolder,'FalseAlarmRate',0.1,'NumCascadeStages',5);
toc
%%
% Use the newly trained classifier to detect a stop sign in an image.
detector = vision.CascadeObjectDetector('CascadePeople.xml');
%%
% Read the test image.
img = imread('~/Documents/MATLAB/Pre/data/SMSQ20/frame-1.jpg');
%%
% Detect a stop sign.
bbox = step(detector,img); 
%%
% Insert bounding box rectangles and return the marked image.
 detectedImg = insertObjectAnnotation(img,'rectangle',bbox,'person');
%%
% Display the detected stop sign.
figure; imshow(detectedImg);
%%
% Remove the image directory from the path.
rmpath(imDir); 

load handel
sound(y,Fs);