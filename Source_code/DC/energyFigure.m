%%
% run several randruns for energy plot

global LOG_allens
allscen=[42 23 25 27 71 72];
% allscen=[23 25];
allinf=[1 2 4 7 8 10];
% allinf=[8 9 10];
allrandruns=1:5;

infParams(1).infParam=[100,0,0]; % TRWS: nIter, randStart, doBPS
infParams(2).infParam=[1,1,0]; % QPBO: useProb, strongPer, useImpr
infParams(3).infParam=[500,1,1,1]; % MQPBO: rounds, stronPers, Kovt., perm.
infParams(4).infParam=[]; % ICM
infParams(5).infParam=[500,1,0]; % TRWSi: maxIter, fastComp, absPrec
infParams(6).infParam=1000; % FastPD: numIter
infParams(7).infParam=[100,1e-7,0]; % BP: numIter, convBound, damping
infParams(8).infParam=[100,1e-7,0]; % TRBP: numIter, convBound, damping
infParams(9).infParam=[0,2]; % LazyFlipper: multilabel, maxSubgraph
infParams(10).infParam=[100,0,1]; % TRWS: nIter, randStart, doBPS
infAlgNames={'TRWS','QPBO','MQPBO','ICM','TRWSi','FastPD','BP','TRBP','LF','LBP'};

opt=getDCOptions;
opt.frames=1:50;
opt.maxItCount=50;



%% discrete
% fix continuous
opt.conOpt.alg='CGD'; % CGD, fmincon
opt.conOpt.maxIter=200;
opt.conOpt.initSimple=1;

for s=allscen
    logsstruct=[];
    scenario=s;
    for ia=allinf
        for r=allrandruns
            disp([s ia r]);
            opt.disOpt.alg=ia;
            if ia==10, opt.disOpt.alg=1; end % LoopyBP
            opt.disOpt.infParam=infParams(ia).infParam;
            opt.randrun=r;
            [metrics2d, metrics3d, allens, stateInfo]=dcTracker(scenario,opt);
            fname=sprintf('data/s%d-r%d-disalg%d-conalg%s.mat', ...
                scenario,opt.randrun,ia,opt.conOpt.alg);
            save(fname,'*allens','metrics*','stateInfo');
            logsstruct(ia,r).log=LOG_allens;
        end
    end
    save(sprintf('d:/visinf/papers/2013/milan-pami/data/enlogsDiscrete-s%d.mat',scenario),'logsstruct');
end


%% continuous
allinf=1:8;
infAlgNames={'CGD','Trust-region','L-BFGS','Simplex','CGD*','Trust-region*','L-BFGS*','Simplex*'};

%fix discrete
opt.disOpt.alg=1; % 1=TRWS, 2=QPBO, 3=MQPBO, 4=ICM, 5=TRWSi, 6=FastPD
opt.disOpt.infParam=[100,0,0]; % TRWS: nIter, randStart, doBPS

% allscen=23;
% opt.frames=1:30;
allscen=[25];
for s=allscen
    scenario=s;
    logsstruct=[];
    for ia=allinf
        switch (ia)
            case {1,5}
                opt.conOpt.alg='CGD'; % CGD, fmincon
                opt.conOpt.maxIter=200;
            case {2,6}
                opt.conOpt.alg='fminunc'; % CGD, fmincon
                opt.conOpt.maxIter=5;
                opt.conOpt.optimoptions= ...
                    optimoptions('fminunc','MaxIter',opt.conOpt.maxIter,'MaxFunEvals',100,'GradObj','on','Display','off');
            case {3,7}
                opt.conOpt.alg='LBFGS'; %
                opt.conOpt.maxIter=400;
                opt.conOpt.optimoptions= ...
                    struct('GradObj','on','Display','off','MaxIter',opt.conOpt.maxIter,'MaxFunEvals',100, ...
                    'LargeScale','off','HessUpdate','bfgs','InitialHessType','identity','GoalsExactAchieve',1,'GradConstr',false);
            case {4,8}
                opt.conOpt.alg='simplex'; % Simplex scheme
                opt.conOpt.maxIter=400;
                opt.conOpt.optimoptions= ...
                    optimset('MaxIter',opt.conOpt.maxIter,'MaxFunEvals',100,'Display','off');
        end
        if ia>4
            opt.conOpt.initSimple=1;
        else
            opt.conOpt.initSimple=0;
        end
        
        for r=allrandruns
            disp([s ia r]);
            opt.randrun=r;
            fname=sprintf('data/s%d-r%d-disalg%d-conalg%s-is%d.mat', ...
                scenario,opt.randrun,opt.disOpt.alg,opt.conOpt.alg,opt.conOpt.initSimple);
            if exist(fname,'file')
                load(fname);
            else
                [metrics2d, metrics3d, allens, stateInfo]=dcTracker(scenario,opt);
                save(fname,'*allens','metrics*','stateInfo');                
            end
            logsstruct(ia,r).log=LOG_allens;
        end
    end
    save(sprintf('d:/visinf/papers/2013/milan-pami/data/enlogsContinuous-s%d.mat',scenario),'logsstruct');
end
