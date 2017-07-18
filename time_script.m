%%% Time Script %%%
clear workspace

%% Handling time
str = date;
fileID = fopen('timeData.txt', 'at');
fprintf(fileID, "Test %s \n", str);

%% Choose programm

var_MOT = false;
var_MOT_v2 = false;
var_detect_people = true;
var_MOT_KLT = false;
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
    filename = 'data/SMSQ20.m4v';
end

%% Displaying results

display = true;

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
    toc
    fprintf(fileID, "algo = %s , display = %d, on video = %s, execution Time = %d \n", "KLT", display, filename, execTime);
end
fclose(fileID);