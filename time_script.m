%%% Time Script %%%
clear workspace

%% Handling time
str = date;
fileID = fopen('timeData.txt', 'at');
fprintf(fileID, "Test %s \n", str);

%% Choose programm

var_MOT = false;
var_MOT_v2 = false;
var_detect_people = false;
var_MOT_KLT = true;
%% Choose video

var_MOT16 = false;
var_atrium = false;
var_Smsq1 = false;
var_Smsq2 = false;
var_SMSQ13 = false;
var_SMSQ20 = true;

if var_MOT16
    filename = 'data/MOT16-04.mp4';
end
if var_atrium
    filename = 'atrium.mp4';
end
if var_Smsq1
    filename = 'data/SmSq1.m4v';
end
if var_Smsq2
    filename = 'data/SmSq2.m4v';
end
if var_SMSQ13
    filename = 'data/SMSQ13.MOV';
end
if var_SMSQ20
    filename = 'data/SMSQ20short50.avi';
end

%% Displaying results

display = false;

%% Let's go!


if var_MOT
    tic;
    memory = MOT(filename, display);
    execTime = toc;
    fprintf(fileID, "algo = %s , display = %d, on video = %s, execution Time = %d \n", "MOT", display, filename, execTime);
end
if var_MOT_v2
    tic;
    memory2 = MOT_v2(filename, display);
    execTime = toc;
    fprintf(fileID, "algo = %s , display = %d, on video = %s, execution Time = %d \n", "MOT v2", display, filename, execTime);
end
if var_detect_people
    tic;
    detect_people(filename, display);
    fprintf(fileID, "algo = %s , display = %d, on video = %s, execution Time = %d \n", "detection", display, filename, execTime);
end
if var_MOT_KLT
    tic;
    memory3 = MOT_KLT(filename, display);
    execTime = toc;
    end1 = 1;
    for i = 1 : size(memory3,2) % id
        l = size(memory3(i).bboxes,1);
        fr(end1:end1+l-1,2)= i * ones(l,1);
        fr(end1:end1+l-1,1)= memory3(i).frame;
        fr(end1:end1+l-1,3:6)= memory3(i).bboxes;
        fr(end1:end1+l-1,7)= 1 * ones(l,1);
        
        end1 = end1 + l;
    end
    csvwrite('resultsKLT.txt', fr);
    
    
    fprintf(fileID, "algo = %s , display = %d, on video = %s, execution Time = %d \n", "KLT", display, filename, execTime);
end
fclose(fileID);