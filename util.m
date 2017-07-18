% Convert detections from txt files to table for trainingImageLabeler



% M = csvread('/Users/CamilleT/Documents/MATLAB/Pre/Source_code/ACF/results/SMSQ20short.txt');
% Dir = dir('/Users/CamilleT/Documents/MATLAB/Pre/data/SMSQ20/');
% inds = [Dir(:).isdir]==0;
% Dir = Dir(inds); 
% [~, index] = natsortfiles({Dir.name}) ;
% Dir = Dir(index);
% BB = cell(1,M(end,1));
% 
% for i = 1 : M(end,1)
%     BB{i} = M(M(:,1)==i,3:6);
% end
% Table1 = table(BB');
% Table2 = struct2table(Dir);
% 
% Table3 = Table2(2:200,1);
% 
% Table3 = [Table3,Table1];

% Convert trainingImageLabeler output to txt file

% %TableDet = output from trainingImageLabeler
% end1 = 503;
% for i = 26 : height(TableDet)
%     l = size(cell2mat(TableDet.Unnamed(i)),1);
%     bb(end1:end1+l-1,1) = i * ones(l,1);
%     bb(end1:end1+l-1,2) = -1 * ones(l,1);
%     bb(end1:end1+l-1,3:6) = cell2mat(TableDet.Unnamed(i));
%     end1 = end1 + l;
% end
    
    

% % Display one by one the detection in one frame
% I = imread(fullfile(cd, 'data', 'SMSQ20','frame-49.jpg'));
% for i = 1015 : 1036
%     Idisp = insertObjectAnnotation(I, 'rectangle', bb(i,3:6), i);
%     imshow(Idisp);
%     pause;
% end
% 
% % Change ids in the matrix of the ground truth
% 
% bb(993 : 1014, 2) = [10 19 13 18 5 21 22 9 1 2 3 14 4 6 20 12 8 7 17 11 16 15]';
%   
%  %Check if everything is well labeled
%  for i = 40:49
%      i
%     I = imread(fullfile(cd, 'data', 'SMSQ20',strcat('frame-', num2str(i), '.jpg')));
%     bbox = bb(bb(:,1) == i+1, 3:6);
%     lab = bb(bb(:,1) == i+1, 2);
%     Idisp = insertObjectAnnotation(I, 'rectangle', bbox, lab);
%     imshow(Idisp);
%     pause;
%  end


% Change the number of frames of a video

% video = VideoReader('/Users/CamilleT/Documents/MATLAB/Pre/data/SMSQ20.m4v');
% outputVideo = VideoWriter('/Users/CamilleT/Documents/MATLAB/Pre/data/SMSQ20short50.avi');
% open(outputVideo);
% 
% for ii = 1:50
%     frame = readFrame(video);
%     writeVideo(outputVideo,frame); 
% end
% 
% close(outputVideo);
