% Convert detections from txt files to table for trainingImageLabeler


% 
% M = csvread('/Users/CamilleT/Documents/MATLAB/Pre/Source_code/ACF/results/SMSQ17.txt');
% Dir = dir('/Users/CamilleT/Documents/MATLAB/Pre/data/SMSQ17/images');
% inds = [Dir(:).isdir]==0;
% Dir = Dir(inds); 
% [~, index] = natsortfiles({Dir.name}) ;
% Dir = Dir(index);
% BB = cell(1,50);
% 
% for i = 1 : 50
%     BB{i} = M(M(:,1)==i,3:6);
% end
% Table1 = table(BB');
% Table2 = struct2table(Dir);
% Table2.name = strcat(Table2.folder, '/', Table2.name); 
% Table3 = Table2(2:51,1);
% 
% Table3 = [Table3,Table1];

% Convert trainingImageLabeler output to txt file

%TableDet = output from trainingImageLabeler
% end1 = 1;
% for i = 1 : height(TableDet)
%     i
%     X = cell2mat(TableDet.Unnamed(i));
%     l = length(X);
%     bb(end1:end1+l-1,1) = i * ones(l,1);
%     bb(end1:end1+l-1,2) = -1 * ones(l,1);
%     bb(end1:end1+l-1,3:6) = cell2mat(TableDet.Unnamed(i));
%     end1 = end1 + l;
% end
%     
%     
% 
% % Display one by one the detection in one frame
% I = imread(fullfile(cd, 'data', 'SMSQ17','images','0051.jpg'));
% for i = 181 : 200
%     Idisp = insertObjectAnnotation(I, 'rectangle', bb(i,3:6), i);
%     imshow(Idisp);
%     pause;
% end
% close figure 1

% 
% % Change ids in the matrix of the ground truth
% 
% bb(1: 4, 2) = [1 2 3 4]';
%   
%  %Check if everything is well labeled
%  for i = 10:15
%      i
%     I = imread(fullfile(cd, 'data', 'SMSQ17','images',strcat('00', num2str(i), '.jpg')));
%     bbox = bb(bb(:,1) == i+1, 3:6);
%     lab = bb(bb(:,1) == i+1, 2);
%     Idisp = insertObjectAnnotation(I, 'rectangle', bbox, lab);
%     imshow(Idisp);
%     pause;
%  end

%close figure 1
% Change the number of frames of a video
% 
video = VideoReader('/Users/CamilleT/Documents/MATLAB/Pre/data/SMSQ17.MOV');
outputVideo = VideoWriter('/Users/CamilleT/Documents/MATLAB/Pre/data/SMSQ17short.avi');
open(outputVideo);
for ii = 1 : 5
    readFrame(video);
end
for ii = 6:55
    frame = readFrame(video);
    writeVideo(outputVideo,frame); 
end

close(outputVideo);
