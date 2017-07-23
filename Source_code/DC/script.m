%%% Script to lauch code 2 %%%

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

csvwrite('result.txt', fr);
% fr = frame number , id , x , y , w , h


% % Try to save the video of the results
% workingDir = '~/Documents/MATLAB/Pre/data/SMSQ13';
% imageNames = dir(fullfile(workingDir,'images','*.jpg'));
% imageNames = {imageNames.name}';
% outputVideo = VideoWriter(fullfile(workingDir,'SMSQ13results.avi'));
% open(outputVideo)
% for ii = 1 : numel(imageNames)
%     img = imread(fullfile(workingDir,'images',cell2mat(imageNames(ii)))); 
%     %fullfile(workingDir,'images',cell2mat(imageNames(ii)))
%     ind = fr(:,1) == ii;
%     bb = fr(ind,3:6);
%     lab = fr(ind,2);
%     inds = [bb(:,3)]>0;
%     bb = bb(inds, :);
%     lab = lab(inds);
%     if ~isempty(bb)
%         for i= 1 : size(bb,1)
%         img = insertObjectAnnotation(img, 'rectangle', bb(i, :), lab(i));
%         end
%         imshow(img);
%     end
%     writeVideo(outputVideo,img)
%     ii
% end
% 
% close(outputVideo)

% Display
displayTrackingResult(stateInfo.sceneInfo, stateInfo);
