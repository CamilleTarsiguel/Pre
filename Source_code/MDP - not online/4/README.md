# The Way They Move: Tracking Targets with Similar Appearance
*Caglayan Dicle, Octavia Camps, Mario Sznaier*


This is the official code repository for our [ICCV 2013 paper](http://www.cv-foundation.org/openaccess/content_iccv_2013/papers/Dicle_The_Way_They_2013_ICCV_paper.pdf)

    @inproceedings{Dicle2013iccv,
     author = {Dicle, Caglayan and Camps, Octavia and Sznaier, Mario},
     title = {The Way They Move: Tracking Targets with Similar Appearance},
     booktitle = {ICCV},
     year = {2013},
    }

### How to run the SMOT code with SMOT data?

 (1) Download the code. It will create a smot directory.
 
 (2) Download the data from [here](http://coe.neu.edu/~cdicle/data/smot_data.zip) and extract it in smot directory. so that you have directory structure like,
 	
 	smot>
 		smot_core
 		smot_data
 		smot_test
 		smot_util
 		 	
 (3) Check if the dataset path is set correctly in file smot/smot_test/test_smot_batch.m
 
 (4) Run test_smot_batch.m

 (5) Please cite to the paper
 	
	
 	
### Prerequisites
The code works on MATLAB R2012a and that is pretty much the only requirement. I did not test the code on other versions, however I expect it to work on newer versions, too, since I did not use obsolete functionality.

Optionaly, if you want to test Interior Points (IP) method and have "patience" to run experiments with IP method, you need to install CVX toolbox from [CVX website](http://www.cvxr.com)


### How to run the SMOT code on your data?

This code merges the points that are input to the algorithm in idl structure. 

idl is an struct array with field xy which holds the detection centers for each frame.

		idl(t).xy   : Nx2 center xy locations  

Here N is number of detections at time t. 

As one can expect parameters are very important to get good results. param is also structure with following fields

		param.similarity_method : rank estimation method 
		min_s : minimum similarity threshold to merge two tracklets
		hor   : horizon window that the algorithm operate in batches. Say if this is 60. Algorithm will process frames in batches 1-60, 61-120, 121-180 …
		eta_max: inlier noise estimate for the detector. 
		debug : debug flag. just set 0. 
		
The most crucial parameters are eta_max and horizon. One needs to play with them to get good results. Please see initialize_smot.m file to get an idea of how to set those parameters. 

If you have the idl and param ready, you can call the algorithm then by typing

		[itl,etime] = smot_associate(idl,param);

The algorithm outputs itl structure with following fields 
		
		itl(n).t_start : start time of track n
		itl(n).t_end   : end time of track n
		itl(n).length  : length of track n
		itl(n).omega   : 1/0 indicator array for occlusion
		itl(n).data    : 2xlength xy locations of track n
				
where n is the index of the track. You can treat it as the id of the track,too. 

see smot_test/test_smot_demo for a simple example.

### License
MIT License except LAPJV functions see LICENSE file.