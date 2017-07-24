function stateInfo=runDP(detections, opt, myopt)

global sceneInfo

% check if opt struct is what we want
if ~isfield(opt,'c_en')
    opt=getPirOptions(opt);
end

stateInfo=getStateForEmptySolution(sceneInfo,opt);
F = length(sceneInfo.frameNums);
frames=1:F;

% frames=1:length(sceneInfo.frameNums);
% % frames=1:100;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % frames=1:30; % do a part of the whole sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if isfield(myopt,'frames'), frames=myopt.frames; end
% 
% if length(sceneInfo.frameNums)<frames(end), frames=1:length(sceneInfo.frameNums); end % if section outside range
% sceneInfo.frameNums=sceneInfo.frameNums(frames);

% myopt.seqLength=length(frames);

dres=detectionsToPirsiavash(detections);

% datadir  = 'data/';
cachedir = 'cache/';
if ~exist(cachedir,'dir'), mkdir(cachedir);end
vid_name = sceneInfo.sequence;
% vid_path = sceneInfo.imgFolder;


%%% Adding transition links to the graph by fiding overlapping detections in consequent frames.
display('in building the graph...')
% fname = [cachedir vid_name '_graph_res.mat'];
dres = build_graph(dres);

% try
%   load(fname)
% catch
%   dres = build_graph(dres);
%   save (fname, 'dres');
% end




%%% setting parameters for tracking
c_en      = .1;     %% birth cost
c_ex      = .1;     %% death cost
c_ij      = .1;      %% transition cost
betta     = 0.2;    %% betta
max_it    = 250;    %% max number of iterations (max number of tracks)
thr_cost  = 18;     %% max acceptable cost for a track (increase it to have more tracks.)


% defaults
% c_en      = 5;     %% birth cost
% c_ex      = 5;     %% death cost
% c_ij      = .1;      %% transition cost
% betta     = 0.2;    %% betta
% max_it    = Inf;    %% max number of iterations (max number of tracks)
% thr_cost  = 18;     %% max acceptable cost for a track (increase it to have more tracks.)

if nargin>1
    c_en=opt.c_en;
    c_ex=opt.c_ex;
    c_ij=opt.c_ij;
    betta=opt.betta;
    thr_cost=opt.thr_cost;
end


display('in DP tracking ...')
tic
% dres_dp       = tracking_dp(dres, c_en, c_ex, c_ij, betta, thr_cost, max_it, 0);
% dres_dp.r     = -dres_dp.id;
toc

tic
display('in DP tracking with nms in the loop...')
dres_dp_nms   = tracking_dp(dres, c_en, c_ex, c_ij, betta, thr_cost, max_it, 1);
dres_dp_nms.r = -dres_dp_nms.id;
toc

% tic
% display('in push relabel algorithm ...')
% dres_push_relabel   = tracking_push_relabel(dres, c_en, c_ex, c_ij, betta, max_it);
% dres_push_relabel.r = -dres_push_relabel.id;
% toc
%%% We ignore the first frame in evaluation since there is no ground truth for it.
% dres              = sub(dres,               find(dres.fr              >1));
% dres_dp           = sub(dres_dp,            find(dres_dp.fr           >1));
% dres_dp_nms       = sub(dres_dp_nms,        find(dres_dp_nms.fr       >1));
% dres_push_relabel = sub(dres_push_relabel,  find(dres_push_relabel.fr >1));

%%
fnum = max(dres.fr);
bboxes_tracked = dres2bboxes(dres_dp_nms, fnum);  %% we are visualizing the "DP with NMS in the lop" results. Can be changed to show the results of DP or push relabel algorithm.
% quick hack
% if scenario==72
%     bboxes_tracked(201).bbox=[];
% end
%% pad rest
if length(bboxes_tracked)<F
    for pp=length(bboxes_tracked)+1:F
        bboxes_tracked(pp).bbox=[];
    end
end



stateInfo=boxesToStateInfo(bboxes_tracked,sceneInfo);

stateInfo=postProcessState(stateInfo);
stateInfo.sceneInfo=sceneInfo;
stateInfo.opt = myopt;
% myopt.conOpt
if isfield(myopt,'dataFunction')
  myopt=getAuxOpt('aa',myopt,sceneInfo,stateInfo.F);
  alldpoints=createAllDetPoints(detections);
  stateInfo.splines=getSplinesFromGT(stateInfo.X,stateInfo.Y,frames,alldpoints,stateInfo.F);
  stateInfo.labeling=getGreedyLabeling(alldpoints,stateInfo.splines,stateInfo.F);
  stateInfo.outlierLabel=length(stateInfo.splines)+1;
  stateInfo.opt = myopt;
end

