Discrete-Continuous Multi-target Tracking
===========================================

This is a framework for multiple target tracking by discrete-continuous
energy minimization. The main idea was first described in this CVPR 2012 paper [(pdf)](http://www.milanton.de/files/cvpr2012/cvpr2012-anton.pdf)


    Discrete-Continuous Optimization for Multi-Target Tracking
    A. Andriyenko, K. Schindler and S. Roth, In: CVPR 2012    

and later extended to include exclusion and statistically derived energy potentials [(pdf)](http://www.milanton.de/files/cvpr2013/cvpr2013-anton.pdf):

    Detection- and Trajectory-Level Exclusion in Multiple Object Tracking
    A. Milan, K. Schindler and S. Roth, In: CVPR 2013    
	
A consolidated manuscript has been submitted to PAMI:

    Multi-Target Tracking by Discrete-Continuous Energy Minimization
    A. Milan, K. Schindler and S. Roth. Submitted to IEEE PAMI.


Installation
==============

This section describes how to get dctracking running.

Get the latest version of the code and cd into that directory

    hg clone https://bitbucket.org/amilan/dctracking
    cd dctracking
    
### OpenGM2 (Linux only)
You can skip this step if you only want to run our simplified model from CVPR 2012. To that end, you can simply remove to energy weights: exclusionFactor=0 and proxcostFactor=0 in the options file.

However, for inference with non-submodular energies, you will need to install the OpenGM2 framework.
In particular, you should download and build the code with QPBO and/or TRW-S support.

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/external/opengm/BUILD/src/external/

    
    
        
Running
=======

Start Matlab and run

	installDCTracker;
	
Now all should be set up. You can start the tracker with.

    [metrics2d, metrics3d, allen, stateInfo, sceneInfo]=swDCTracker('scenes/TUD-Campus.ini','config/default2dSimple.ini');
    
The output should be similar to the one in logs/log_tud-campus.txt. Note that both discrete inference 
and continuous minimization may lead to slightly different final results, depending on the current 
software environment and hardware architecture.
    
You can display the results by calling

    displayTrackingResult(stateInfo.sceneInfo, stateInfo);
    
    
Other videos
------------

To run the tracker on other videos, adjust the necessary settings in a scene file in 

    ./scenes/... .ini
    
Parameters can be set in

    ./config/... .ini

	
Please do not forget to cite our work if you end up using this code:

    @inproceedings{Milan:2013:DTE,
	    Author = {Anton Milan and Konrad Schindler and Stefan Roth},
	    Booktitle = {CVPR},
	    Title = {Detection- and Trajectory-Level Exclusion in Multiple Object Tracking},
	    Year = {2013}
    }

	@inproceedings{Andriyenko:2012:DCO,
		Author = {Anton Andriyenko and Konrad Schindler and Stefan Roth},
		Booktitle = {CVPR},
		Title = {Discrete-Continuous Optimization for Multi-Target Tracking},
		Year = {2012}
	}	