function detect_people(filename,videoname, detector, display)
 % Use various detectors to detect pedestrians
 % output : a txt file with the detections as follow
 %frame; -1; x; y; w; h; confidence score; (ignored/considered; type of
 %object; (visibility ratio))
root = cd;
%% Choose Detector
switch detector
    case 'ACF'
    var_peopleDetectorACF = true;
    var_CascadeObjectDetector = false;
    var_PeopleDetector = false;
    var_DPMDetector = false;
    case 'Cascade'
    var_peopleDetectorACF = false;
    var_CascadeObjectDetector = true;
    var_PeopleDetector = false;
    var_DPMDetector = false;
    case 'HOG'
    var_peopleDetectorACF = false;
    var_CascadeObjectDetector = false;
    var_PeopleDetector = true;
    var_DPMDetector = false;
    case 'DPM'
    var_peopleDetectorACF = false;
    var_CascadeObjectDetector = false;
    var_PeopleDetector = false;
    var_DPMDetector = true;
    otherwise
        error('This detector is not supported');
end

if var_peopleDetectorACF
    str = 'ACF';
end

if var_PeopleDetector
    str = 'HOG';
end

if var_CascadeObjectDetector
    str = 'Cascade';
end
if var_DPMDetector
    str = 'DPM';
end

str = fullfile(root,'Source_code', str, 'results', strcat(videoname ,'.txt'));
memory = [];

%% Choose parameters

maxH = 100;
minH = 2;
maxW = 70;
minW = 2;
%% Create System objects used for reading video, detecting moving objects,
% and displaying the results.
obj = setupSystemObjects();


% Go to the first frame where there is a detection

frame = readframe();
numframe = 1
[score, bboxes] = detectObjects(frame);
while isempty(bboxes)
    frame = readframe();
    numframe = numframe + 1
[score, bboxes] = detectObjects(frame);
end

% Detect moving objects, and track them across video frames.
while hasFrame(obj.reader)
    
    [score, bboxes] = detectObjects(frame);
    ind = size(memory, 1);
    memory(ind + 1 : ind + length(score), 1) = numframe * ones(length(score),1);
    memory(ind + 1 : ind + length(score), 2) = -1 * ones(length(score),1);
    memory(ind + 1: ind + length(score), 3:6) = bboxes;
    memory(ind + 1: ind + length(score), 7) = score;
    
    if display
        displayDetectingResults();
    end
   frame = readframe();
   numframe = numframe + 1
end

csvwrite(str, memory);
%% Create System Objects
% Create System objects used for reading the video frames, detecting
% foreground objects, and displaying results.

    function obj = setupSystemObjects()
        % Initialize Video I/O
        % Create objects for reading a video from a file, drawing the tracked
        % objects in each frame, and playing the video.
        
        % Create a video file reader.
        obj.reader = VideoReader(filename);
        
        % Create two video players, one to display the video,
        % and one to display the foreground mask.
        obj.videoPlayer = vision.VideoPlayer('Position', [20, 400, 700, 400]);
        
        % Create System objects for foreground detection and blob analysis
        
       
        if var_peopleDetectorACF
            obj.detector = peopleDetectorACF('caltech-50x21');
        end
        if var_CascadeObjectDetector
            obj.detector = vision.CascadeObjectDetector('UpperBody');
        end
        if var_PeopleDetector
            obj.detector = vision.PeopleDetector('UprightPeople_96x48');
        end
        if var_DPMDetector
            load(fullfile(cd,'Source_code','DPM','PROJECT/2007/car_final'));
            model.vis = @() visualizemodel(model, ...
                  1:2:length(model.rules{model.start}));
            obj.detector = model;
        end
        
    end


%% Read a Video Frame
% Read the next video frame from the video file.
    function frame = readframe()
        frame = readFrame(obj.reader);
    end

%% Detect Objects
% The |detectObjects| function returns the centroids and the bounding boxes
% of the detected objects. It also returns the binary mask, which has the 
% same size as the input frame. Pixels with a value of 1 correspond to the
% foreground, and pixels with a value of 0 correspond to the background.   
%
% The function performs motion segmentation using the foreground detector. 
% It then performs morphological operations on the resulting binary mask to
% remove noisy pixels and to fill the holes in the remaining blobs.  

    function [score, bboxes] = detectObjects(frame)

        if var_peopleDetectorACF
            [bboxes, score ] = detect(obj.detector, frame);
        end
        if var_CascadeObjectDetector
            bboxes = step(obj.detector, frame);
            score = ones(size(bboxes,1),1);
        end
        
        if var_PeopleDetector
            [bboxes, score ] = step(obj.detector, frame);
        end
        if var_DPMDetector
            [bboxes,score] = test(frame, obj.detector, 100);
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
    end




%% Display Detecting Results
% The |displayTrackingResults| function draws a bounding box and label ID 
% for each track on the video frame and the foreground mask. It then 
% displays the frame and the mask in their respective video players. 

    function displayDetectingResults()
        % Convert the frame and the mask to uint8 RGB.
        %frame = im2uint8(frame);

                frame = insertObjectAnnotation(frame, 'rectangle', ...
                  bboxes, score);
        
        % Display the frame.
          
        obj.videoPlayer.step(frame);
    end


%% Test function for DPM
function [bboxes,score] = test(frame, model, num_dets)
cls = model.class;
%fprintf('///// Running demo for %s /////\n\n', cls);

% detect objects

[ds, bs] = imgdetect(frame, model, -1);
top = nms(ds, 0.5);
top = top(1:min(length(top), num_dets));
ds = ds(top, :);
bs = bs(top, :);

if model.type == model_types.Grammar
  bboxes = [ds(:,1:4) bs];
end
%showboxes(frame, reduceboxes(model, bs)); % I don't care about that


if model.type == model_types.MixStar
  % get bounding boxes
  bbox = bboxpred_get(model.bboxpred, ds, reduceboxes(model, bs));
  bbox = clipboxes(frame, bbox);
  top = nms(bbox, 0.5);
  bboxes = bbox(top,:);
end
score = ones(size(bboxes,1),1);
end

end



