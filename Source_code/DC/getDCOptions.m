function opt=getDCOptions()
% create a struct with all necessary options and parameters

global scenario;


% tracking on image or on ground plane?
opt.track3d=howToTrack(scenario);

% general
opt.verbosity=3;                % 0=silent, 1=short info, 2=long info, 3=all
opt.visOptim=0;                 % visualize optimization
opt.met2d=1;                    % always compute metrics in 2d (slower)
opt.maxModels=20000;            % max number of trajectory hypotheses
opt.keepHistory=2;              % keep unused models for n iterations
opt.cutToTA=0;                  % cut detections, ground truth and result to tracking area
opt.randrun=1;                  % random seed
opt.remOcc=0;                   % remove occluded GT and result
opt.maxItCount=100;             % abort after max iterations reached
opt.occ=0;                      % turn on / off occlusion reasoning

% unary (data) term, 
opt.dataFunction = 2;           % 1 = L2-dist, 2=L2-dist^2, 3=dist Lorenzian

% curvature
opt.speedFunction = 2;          % 1 = max(...),  2 = curvature, 3 = curv. lor.

% fidelity
opt.fidFunction = 3;            % 1 = original (^3), 2 = linear, 3 = lorenzian

opt.nInitModels= 10;			% initial hypotheses (Sec. 5.1)

% models modifications (Sec. 5.2)
opt.nMaxAddExt=50;
opt.nMaxAddMerged=50;
opt.nAddRandomModels=100;
opt.nAddModelsFromOutliers=100;

% all equal
newmodels=10;
opt.nMaxAddExt=newmodels;opt.nMaxAddMerged=newmodels;opt.nAddRandomModels=newmodels;opt.nAddModelsFromOutliers=newmodels;


opt.minCPs=         1;
opt.ncpsPerFrame=   1/10;
opt.tau =           10;     % threshold (pixel) for spatio-temporal neighbors
opt.borderMargin =  100;   % (pixel) % distance for persistence

opt.startFromEKF=1:5;
opt.startFromPir=1;
opt.startFromGT=0;
opt.EKFDir='data/init/ekftracking';
opt.DPDir='data/init/dptracking';

% Discard all detections below threshold
opt.detThreshold=.0;

% how do we scale detections?
opt.detScale.sigA=0; % shift
opt.detScale.sigB=1; % peakiness



opt.outlierCost=    300;
opt.labelCost =     20;
opt.goodnessFactor= 1;
opt.unaryFactor=   10;



    opt.mex = 1;
    opt.segFactor=1000;
    opt.readjustSE=0;
    
    % PAMI
    opt.labelCost=200;
    opt.persistenceFactor=0;
    opt.dataFunction = 4;           % 1 = L1 Norm, 2=L2 Norm, 3=Lorenzian
    opt.speedFunction = 3;          % 1 = max(...),  2 = curvature, 3 = curv. lor.
    opt.fidFunction = 4;            % 1 = original (^3), 2 = linear, 3 = lorenzian
    opt.curvatureFactor=0;
    opt.slopeFactor=1e-3;
    opt.proxcostFactor=1000;
	opt.exclusionFactor=200;
    opt.fidelityFactor=0;
    opt.pairwiseFactor=1;
    
    opt.frames=1:2000;    
    opt.visOptim=0;
    
    % KITTI
    opt.persistenceFactor=0;
    opt.labelCost=800;
    opt.outlierCost=    800;
    opt.segFactor=10;
    opt.slopeFactor=1e-4;
%     opt.visOptim=1;

if opt.track3d
    opt.tau =           750;   % threshold (mm) for spatio-temporal neighbors
    opt.borderMargin = 2000;    % distance for persistence

    % ogm
    opt.visOptim=0;
    opt.outlierCost=    200;
    opt.labelCost =     1;
    opt.unaryFactor=    100;
    opt.pairwiseFactor=1;
    opt.goodnessFactor= 1;        
    opt.fidelityFactor = 1;
    opt.curvatureFactor = 1;    
    opt.persistenceFactor = .1;
    opt.exclusionFactor=200;
    newmodels=10;
    opt.nMaxAddExt=newmodels;opt.nMaxAddMerged=newmodels;opt.nAddRandomModels=newmodels;opt.nAddModelsFromOutliers=newmodels;
    opt.nInitModels=    10;
    opt.tau =           200;
    
    opt.proxcostFactor= 10.0;
    opt.slopeFactor=.00001;


    
    opt.startFromEKF=1:5;
    opt.startFromPir=1;
    opt.startFromCon=0;
    opt.startFromGT=0;
    opt.nInitModels=    50;

    
    opt.frames=141:2000;


    % experimenting with new cont. en.
    opt.outlierCost=    500;
    opt.labelCost=100;
    opt.unaryFactor=100;    
    opt.fidelityFactor = 0;
    opt.slopeFactor = 20;
    opt.curvatureFactor = .1;
    opt.persistenceFactor=20;
    opt.proxcostFactor=100;
    opt.exclusionFactor=100;
    
    % PAMI
    opt.labelCost=0.3026;
    opt.outlierCost=    520.0769;
    opt.unaryFactor=320.3257;    
    opt.pairwiseFactor=0.21;
    opt.goodnessFactor=2.0151;
    opt.fidelityFactor = 0;
    opt.slopeFactor = 13.3286;
    opt.curvatureFactor = 0.0002;
    opt.persistenceFactor=22.1478;
    opt.proxcostFactor=29.6237;
    opt.exclusionFactor=19.4525;
    

    opt.maxItCount=50;

    opt.segFactor=1e-2;
    opt.readjustSE=0;
    
    opt.randrun=2;
    opt.visOptim=1;

    
%     opt.met2d=0;
%     opt.unaryFactor=10;
%     opt.outlierCost=1000;
%     opt.slopeFactor=.005;
%     opt.labelCost=50;
%     opt.persistenceFactor=1;

    

end

    % discrete optimization options
%     opt.disOpt.alg=1; % 1=TRWS, 2=QPBO, 3=MQPBO, 4=ICM, 5=TRWSi, 6=FastPD
%     opt.disOpt.infParam=[100,0,1]; % TRWS: nIter, randStart, doBPS
    
    opt.disOpt.alg='QPBO'; % 1=TRWS, 2=QPBO, 3=MQPBO, 4=ICM, 5=TRWSi, 6=FastPD
    opt.disOpt.infParam=[1,1,1]; % QPBO: useProb, strongPer, useImpr

%     opt.disOpt.alg=3; % 1=TRWS, 2=QPBO, 3=MQPBO, 4=ICM, 5=TRWSi, 6=FastPD
%     opt.disOpt.infParam=[10,1,1,1]; % MQPBO: rounds, stronPers, Kovt., perm.
% 
%     opt.disOpt.alg=4; % 1=TRWS, 2=QPBO, 3=MQPBO, 4=ICM, 5=TRWSi, 6=FastPD
%     opt.disOpt.infParam=[];
%    
%     opt.disOpt.alg=5; % 1=TRWS, 2=QPBO, 3=MQPBO, 4=ICM, 5=TRWSi, 6=FastPD
%     opt.disOpt.infParam=[500,1,0]; % TRWSi: maxIter, fastComp, absPrec

%     opt.disOpt.alg=6; % 1=TRWS, 2=QPBO, 3=MQPBO, 4=ICM, 5=TRWSi, 6=FastPD
%     opt.disOpt.infParam=100; % FastPD: numIter

%      opt.disOpt.alg=7; % BP
%      opt.disOpt.infParam=[100,1e-7,0]; % BP: numIter, convBound, damping

%      opt.disOpt.alg=8; % TRBP
%      opt.disOpt.infParam=[100,1e-7,0]; % TRBP: numIter, convBound, damping

%       opt.disOpt.alg=9; % Lazy Flipper
%       opt.disOpt.infParam=[0,2]; % LF: multilabel, maxSubgraph
    
	
    
    
% continuous minimization options
    opt.conOpt.alg='CGD'; % CGD, fmincon
    opt.conOpt.maxIter=200;
    
%     opt.conOpt.alg='fminunc'; % CGD, fmincon
%     opt.conOpt.maxIter=5;
%     
    opt.conOpt.alg='LBFGS'; % 
    opt.conOpt.maxIter=200;
%  
%     opt.conOpt.alg='lsfit'; % Least squares fit
%
%     opt.conOpt.alg='simplex'; % Simplex scheme
%     opt.conOpt.maxIter=400;


	% start continuous optimization from Least-Squares fit
    opt.conOpt.initSimple=1;
    

	% old or debug parameters, please ignore
    opt.saveEnSS=0; % save energy snapshots
	opt.oriKappa=0;
	opt.oriWeight=0;
	opt.startFromCon=0;
	opt.startFromLinked=0;
    opt.disOpt.initSimple=0;
    opt.readjustSE=0;
    opt.lengthFactor = 0;
	opt.totalSupFactor= 2;
	opt.meanDPFFactor=  2;
	opt.meanDPFeFactor= 2;
    opt.randFit=0;
    
%     %% CVPR 2012
% 
%     % unary (data) term, 
%     opt.dataFunction = 2;           % 1 = L2-dist, 2=L2-dist^2, 3=dist Lorenzian
% 
%     % curvature
%     opt.speedFunction = 1;          % 1 = max(...),  2 = curvature, 3 = curv. lor.
% 
%     % fidelity
%     opt.fidFunction = 1;            % 1 = original (^3), 2 = linear, 3 = lorenzian    
%     opt.conOpt.alg='lsfit'; % Least squares fit
%     opt.readjustSE=0;
%     opt.conOpt.initSimple=0;
%     opt.disOpt.initSimple=0;
%     opt.proxcostFactor=0;
%     opt.exclusionFactor=0;
%     opt.fidelityFactor=1;


end