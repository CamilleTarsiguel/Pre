function detect_people_img(img_dir, display)
 % Use various detectors to detect pedestrians
 % output : a txt file with the detections as follow
 %frame; -1; x; y; w; h; confidence score; (ignored/considered; type of
 %object; (visibility ratio))

%% Choose Detector

var_peopleDetectorACF = true;
var_CascadeObjectDetector = false;
var_PeopleDetector = false;

if var_peopleDetectorACF
    str = 'ACF';
end

if var_PeopleDetector
    str = 'HOG';
end

if var_CascadeObjectDetector
    str = 'Cascade';
end

str = strcat('detections-', str, '.txt');
memory = [];


%% Create System objects used for reading video, detecting moving objects,
% and displaying the results.
obj = setupSystemObjects();


% Go to the first frame where there is a detection
numframe = 1
frame = readFrame();

[score, bboxes] = detectObjects(frame);
while isempty(bboxes)
    numframe = numframe + 1
    frame = readFrame();
    [score, bboxes] = detectObjects(frame);
end

% Detect moving objects, and track them across video frames.
while (numframe < obj.end)
    
    [score, bboxes] = detectObjects(frame);
    ind = size(memory, 1);
    memory(ind + 1 : ind + length(score), 1) = numframe * ones(length(score),1);
    memory(ind + 1 : ind + length(score), 2) = -1 * ones(length(score),1);
    memory(ind + 1: ind + length(score), 3:6) = bboxes;
    memory(ind + 1: ind + length(score), 7) = score;
    
    if display
        displayDetectingResults();
    end
   frame = readFrame();
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
        obj.IMG = dir(img_dir);
        inds = [obj.IMG(:).bytes]>0;
        obj.IMG = obj.IMG(inds); % struct listing all the files in the img_dir folder with info like name , date, bytes, folder, isfolder..
        [~, index] = natsortfiles({obj.IMG.name}) ;
        obj.IMG = obj.IMG(index);
        obj.end = length(obj.IMG);
%         IMG(1).name
       
        if var_peopleDetectorACF
            obj.detector = peopleDetectorACF('caltech-50x21');
        end
        if var_CascadeObjectDetector
            obj.detector = vision.CascadeObjectDetector('UpperBody');
        end
        if var_PeopleDetector
            obj.detector = vision.PeopleDetector('UprightPeople_96x48');
        end
        
    end


%% Read a Video Frame
% Read the next video frame from the video file.
    function frame = readFrame()

         frame = imread(strcat(img_dir, (obj.IMG(numframe).name)));
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
                frame = insertText(frame, [10, 10], num2str(numframe));
        % Display the frame.
          
        imshow(frame);
    end

end



