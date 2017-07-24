%%% Script for loading annotation files %%%

clear files;
clear annotations;
% 
% 
workingDir = '/Users/CamilleT/Documents/MATLAB/Pre/Source_code/dollar/data/Caltech/train04/annotations/';
imageDir = '/Users/CamilleT/Documents/MATLAB/Pre/Source_code/dollar/data/Caltech/train04/images/';
Dir = dir(fullfile(workingDir, '*.txt'));
index = 1;
for i = 1 : numel(Dir)
    filename = Dir(i).name;
    fullname = fullfile(workingDir,filename);
    img_name = fullfile(imageDir, strcat(extractBefore(filename,'.'),'.jpg'));
    data = importdata(fullname);
    if isstruct(data)
        fprintf('Loading annotations for %s...\n',filename);
        files {index} = img_name ;
        annotations{index} = data.data(:,1:4);
        index = index + 1;
    end
end
%trainingTable = table(files', annotations');


