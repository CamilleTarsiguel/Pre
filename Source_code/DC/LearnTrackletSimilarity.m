% try various linking parameters (used for ETH sequences only)
clusternr=4; useexp=1;
maxexper=99;
global gtInfo scenario detections sceneInfo opt
savesim=1;


load([resdir sprintf('res_%03d.mat',useexp-1)]); opts(1).track3d=0; opt=opts(1);
scenario=53;opt.track3d=howToTrack(scenario);sceneInfo=getSceneInfo(scenario);
stateInfo=infos(scenario,bestruns(scenario,useexp)).stateInfo;
saveTrackletSim;

%%
scenario=51;opt.track3d=howToTrack(scenario);sceneInfo=getSceneInfo(scenario);
stateInfo=infos(scenario,bestruns(scenario,useexp)).stateInfo;
saveTrackletSim;
TrackletSimilarityNormFeat;
savesim=0;

%% now compute 
savesim=1;
scenario=53;opt.track3d=howToTrack(scenario);sceneInfo=getSceneInfo(scenario);
stateInfo=infos(scenario,bestruns(scenario,useexp)).stateInfo;

[detections nPoints]=parseDetections(sceneInfo);
% stateInfo=getBBoxesFromState(stateInfo); printFinalEvaluation(stateInfo,gtInfo,sceneInfo,opt);
for exper=0:maxexper, computeTrackletSimilarity; end
    
%
savesim=1;
scenario=51;opt.track3d=howToTrack(scenario);sceneInfo=getSceneInfo(scenario);
stateInfo=infos(scenario,bestruns(scenario,useexp)).stateInfo;
[detections nPoints]=parseDetections(sceneInfo);
% stateInfo=getBBoxesFromState(stateInfo); printFinalEvaluation(stateInfo,gtInfo,sceneInfo,opt);
for exper=0:maxexper, computeTrackletSimilarity; end
savesim=0;

%%
maxexper=99;
% train on bahnhof
fprintf('Training on Bahnhof\n');
scenario=51;
clf;
hold on
allmax=0; allmaxi=0;
for exper=0:maxexper %%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    resfile=sprintf('tmp/linked/simres/%03d/exp%03d-s%02d.mat',clusternr,exper,scenario);
    load(resfile)
    plot(allm2d(:,12),'.-','color',getColorFromID(exper));
    [maxmota atiter]=max(allm2d(:,12));
    if maxmota>allmax
        allmax=maxmota;
        allmaxi=exper;
    end
%     pause
end

bestexp=allmaxi;
resfile=sprintf('tmp/linked/simres/%03d/exp%03d-s%02d.mat',clusternr,bestexp,scenario);
load(resfile)
[maxmota atiter]=max(allm2d(:,12));
upToFactor=allmind(atiter)/firstmin;
fprintf('best exper: %d, theta: %f, upToFactor: %f (at iteration %i)\n',bestexp, allmind(atiter),upToFactor, atiter);

%% test on sunnyday
fprintf('Testing on sunnyday\n');
% useexp=22; resdir='~/visinf/projects/ongoing/dctracking/results/062/';
global gtInfo scenario detections sceneInfo opt
load([resdir sprintf('res_%03d.mat',useexp-1)]); opts(1).track3d=0; opt=opts(1);
scenario=53;opt.track3d=howToTrack(scenario);sceneInfo=getSceneInfo(scenario);
stateInfo=infos(scenario,bestruns(scenario,useexp)).stateInfo;
[detections nPoints]=parseDetections(sceneInfo);
printFinalEvaluation(stateInfo,gtInfo,sceneInfo,opt);
% stateInfo=getBBoxesFromState(stateInfo); printFinalEvaluation(stateInfo,gtInfo,sceneInfo,opt);
exper=bestexp;
% exper=75;
%sim: -0.515678,  scale: 0.998752
%%
computeTrackletSimilarity;
%%
saveYangXML(stInfo,'/home/amilan/diss/others/yangbo/ETH-Person/sunnyday_dcPAMI.xml');
save(sprintf('tmp/linked/results/clnr-%03d-s%d.mat',clusternr,scenario),'stInfo');
evaluateYang;

%%

% train on sunnyday
fprintf('Training on Sunny Day\n');
scenario=53;
clf;
hold on
allmax=0; allmaxi=0;
for exper=0:maxexper
    resfile=sprintf('tmp/linked/simres/%03d/exp%03d-s%02d.mat',clusternr,exper,scenario);
    load(resfile)
%     allmind
    plot(allm2d(:,12),'.-','color',getColorFromID(exper));
    [maxmota atiter]=max(allm2d(:,12));
    if maxmota>allmax
        allmax=maxmota;
        allmaxi=exper;
    end
%     pause
end

bestexp=allmaxi;
resfile=sprintf('tmp/linked/simres/%03d/exp%03d-s%02d.mat',clusternr,bestexp,scenario);
load(resfile)
[maxmota atiter]=max(allm2d(:,12));
upToFactor=allmind(atiter)/firstmin;
fprintf('best exper: %d, theta: %f, upToFactor: %f (at iteration %i)\n',bestexp, allmind(atiter),upToFactor, atiter);
exper=bestexp;

%% test on bahnhof
fprintf('Testing on Bahnhof\n');
% useexp=1; resdir='~/visinf/projects/ongoing/dctracking/results/052/';
% global gtInfo scenario detections sceneInfo opt
load([resdir sprintf('res_%03d.mat',useexp-1)]); opts(1).track3d=0; opt=opts(1);
scenario=51;opt.track3d=howToTrack(scenario);sceneInfo=getSceneInfo(scenario);
stateInfo=infos(scenario,bestruns(scenario,useexp)).stateInfo;
[detections nPoints]=parseDetections(sceneInfo);
printFinalEvaluation(stateInfo,gtInfo,sceneInfo,opt);
% stateInfo=getBBoxesFromState(stateInfo); printFinalEvaluation(stateInfo,gtInfo,sceneInfo,opt);

exper=bestexp;
% exper=75;
%sim: -0.515621,  scale: 0.998831
%%
computeTrackletSimilarity
%%
saveYangXML(stInfo,'/home/amilan/diss/others/yangbo/ETH-Person/bahnhof_dcPAMI.xml');
save(sprintf('tmp/linked/results/clnr-%03d-s%d.mat',clusternr,scenario),'stInfo');

evaluateYang;
!cat ETH-comparison.txt