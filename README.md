# Pre


People Detection and Tracking for Smart Square Project
======================================================

May - July 2017


Author : Camille Tarsiguel, camille.tarsiguel@ensta-paristech.fr

Citation : I used and modified the code from the following articles
- DC Tracking : Discrete-Continuous Optimization for Multi-Target Tracking, A. Andriyenko, K. Schindler and S. Roth, In: CVPR 2012. http://www.milanton.de/dctracking/
- MDP Tracking : Learning to Track: Online Multi-Object Tracking by Decision Making, Xiang, Yu and Alahi, Alexandre and Savarese, Silvio, In: ICCV 2015. https://github.com/yuxng/MDP_Tracking
- ROT Tracking :
- Evaluation : MOT challenge devkit https://motchallenge.net
- ACF detector and Caltech Pedestrian Dataset + Toolbox : P. Dollar
- DPM detector :

This code was tested on Mac OSX 10.12.5 and Matlab R2017a.

INSTALLATION : 
==============

Nothing is needed if working on Mac.
For other types of computer, see individual codes

FOLDER STRUCTURE : 
==================

codes :

main : runs the detectors and trackers —> open file to choose parameters \
main_eval : runs the evaluation

convertVid2img : convert a video to images files (1 image = 1 frame) \
convImg2vid : convert image files from a folder into a video \
detect_people_img : work with the script_detect_people_img \
detect_people : main script for detection (on videos) work with main \
MOT_KLT :  Kanade-Lucas-Tomasi tracker \
MOT_v2 : Kalman based tracker \
MOT : Kalman based tracker with background substraction and blob analysis \
MultiObjectTrackerKLT : class to convert the KLT tracker into an multiple object tracker \
natsort * all these files are helpful to reorder file when calling dir in Matlab \
script_detect_img : run the detect_people_img \
trainCascade : train Cascade detector on Caltech Dataset \
trainRCNN : train RCNN detector on Caltech Dataset \
util : Some short codes that can be useful ( 1. Convert detections from txt files to table for trainingImageLabeler	2. Convert trainingImageLabeler output to txt file	3. Display one by one the detection in one frame and change the ids (for the tracking ground truth)	4. Change the number of frames of a video )



subfolders:

data : # videos from Domplatz for evaluation and display \
	—> SMSQ20  % 50 frames, for evaluation

groundTruth : Contain the groundTruth for evaluation. Each file for gt is a txt file with the name of the video it refers to. Structure is like in MOT Challenge : frame , id , x , y, w , h , (conf)

Source_code : # this is the folder containing the code for the trackers and detectors \
	—> ACF # Agregated channel feature detector (results from MATLAB peopleDetectorACF, no code) \
	—> Cascade # Cascade Detector (results from MATLAB CascadeObjectDetector, no code) \
	—> DC # Discrete continuous tracker (offline) \
	—> dollar # ACF detector \
	—> DPM # Deformable part model detector \
	—> HOG # HOG + SVM detector (results from MATLAB PeopleDetector, no code) \
	—> MDP # Markov Decision Process tracker \
	—> motchallenge-devkit # evaluation \
	—> motutils # some codes for DC \
	—> ROT # Robust online tracker based on tracklet confidence


RUNNING :
=========

To run the code, you only need to open main.m and set the parameters  \
To evaluate (other than running time) open main_eval.m and set the parameters accordingly
