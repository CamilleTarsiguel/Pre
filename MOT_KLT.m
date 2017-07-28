function memory = MOT_KLT(filename, display)


%% Choose detector

var_peopleDetectorACF = true;
var_CascadeObjectDetector = false;
var_PeopleDetector = false;


%% Create System objects used for reading video, detecting moving objects,
% and displaying the results.
obj = setupSystemObjects();



% Go to the first frame where there is a detection
numframe = 1; 
frame = readFrame();

bboxes = detectObjects(frame);
while isempty(bboxes)
    frame = readFrame();
    numframe = numframe + 1;
    bboxes = detectObjects(frame);
    fprintf("Waiting for detection, fr no %d\n", numframe);
end
obj.tracker.addDetections(rgb2gray(frame), bboxes, numframe);
% Detect moving objects, and track them across video frames.
while ~isDone(obj.reader)
    numframe = numframe + 1;
    if mod(numframe, 10) == 0
        bboxes = detectObjects(frame);
        obj.tracker.addDetections(rgb2gray(frame), bboxes, numframe);
    else
        obj.tracker.track(rgb2gray(frame), numframe);
    end
    if display
        displayTrackingResults();
    end
   frame = readFrame();
   
end

 % keep in memory all the tracks
    
memory = obj.tracker.Memory;


%% Create System Objects
% Create System objects used for reading the video frames, detecting
% foreground objects, and displaying results.

    function obj = setupSystemObjects()
        % Initialize Video I/O
        % Create objects for reading a video from a file, drawing the tracked
        % objects in each frame, and playing the video.
        
        % Create a video file reader.
        obj.reader = vision.VideoFileReader(filename);
        
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
        
        obj.tracker = MultiObjectTrackerKLT;
        
    end


%% Read a Video Frame
% Read the next video frame from the video file.
    function frame = readFrame()
        frame = obj.reader.step();
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

    function bboxes = detectObjects(frame)

        maxH = 100;
        minH = 2;
        maxW = 70;
        minW = 2;
        
        if var_peopleDetectorACF
            [bboxes, score ] = detect(obj.detector, frame);
            %centroids = [bboxes(:,1) + 0.5 * bboxes(:,3) bboxes(:,2) - 0.5 * bboxes(:,4)];
        end
        if var_CascadeObjectDetector
            bboxes = step(obj.detector, frame);
            %centroids = [bboxes(:,1) + 0.5 * bboxes(:,3) bboxes(:,2) - 0.5 * bboxes(:,4)];
        end
        
        if var_PeopleDetector
            [bboxes, score ] = step(obj.detector, frame);
            %centroids = [bboxes(:,1) + 0.5 * bboxes(:,3) bboxes(:,2) - 0.5 * bboxes(:,4)];
        end
        
        index = [bboxes(:,3)<maxW];
        bboxes = bboxes(index, :);
       % score = score(index, 1);
        index = [bboxes(:,3)>minW];
        bboxes = bboxes(index, :);
       % score = score(index, 1);
        index = [bboxes(:,4)<maxH];
        bboxes = bboxes(index, :);
       % score = score(index, 1);
        index = [bboxes(:,4)>minH];
        bboxes = bboxes(index, :);
       % score = score(index, 1);
    end


%% Display Tracking Results
% The |displayTrackingResults| function draws a bounding box and label ID 
% for each track on the video frame and the foreground mask. It then 
% displays the frame and the mask in their respective video players. 

    function displayTrackingResults()
        % Convert the frame and the mask to uint8 RGB.
        %frame = im2uint8(frame);
       
        
        displayFrame = insertObjectAnnotation(frame, 'rectangle',...
        obj.tracker.Bboxes, obj.tracker.BoxIds);
        displayFrame = insertMarker(displayFrame, obj.tracker.Points);
        obj.videoPlayer.step(displayFrame);
        
%         
%         
    end

end



