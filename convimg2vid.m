%%% convert image sequence to video

workingDir = 'data/SMSQ13'
imageNames = dir(fullfile(workingDir,'results','*.jpg'));
imageNames = {imageNames.name}';
outputVideo = VideoWriter(fullfile(workingDir,'SMSQ13results.avi'));
outputVideo.FrameRate = shuttleVideo.FrameRate;
open(outputVideo)
for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,'images',imageNames{ii}));
   writeVideo(outputVideo,img)
end
close(outputVideo)

shuttleAvi = VideoReader(fullfile(workingDir,'SMSQ13results.avi'));
ii = 1;
while hasFrame(shuttleAvi)
   mov(ii) = im2frame(readFrame(shuttleAvi));
   ii = ii+1;
end
figure
imshow(mov(1).cdata, 'Border', 'tight')
movie(mov,1,shuttleAvi.FrameRate)