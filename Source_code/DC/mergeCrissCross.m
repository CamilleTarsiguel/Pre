function mhsnew=mergeCrissCross(mhs,alldpoints,T,opt)
[prox, proxt, proxcost, merspl]=getSplineProximity(mhs,T,[],opt);
nmer=size(merspl,1);
addedMerged=0; xi=3;
mhsnew=getEmptyModelStruct();
fitnorm='a';
for nm=1:nmer
    mus=merspl(nm,3);mue=merspl(nm,4);mintime=merspl(nm,5);
    
    mh1=mhs(merspl(nm,1)); mh2=mhs(merspl(nm,2));
    
    tt1=mh1.start:mintime-xi; tt2=mintime+xi:mh2.end;
    
    % fit through splines
    xt1=ppval(mh1,tt1); xt2=ppval(mh2,tt2);
    xy=[xt1 xt2];
    t=[tt1 tt2];
    %             t
    tr=t+0.01*rand(1,length(t))-0.005;
    addedMerged=addedMerged+1;
    sfit=splinefit(tr,xy,round(max(opt.minCPs,(tr(end)-tr(1))*opt.ncpsPerFrame )),fitnorm);
    sfit=adjustSplineStruct(sfit, min(t), max(t), alldpoints, T, 0, [], [], []);
    mhsnew(addedMerged) = sfit;
    
    %             reopenFig('criss cross'); clf; hold on;
    %
    %             allval=ppval(mh1,mh1.start:mh1.end);
    %             allx=allval(1,:);ally=allval(2,:);
    %             plot3(allx,ally,mh1.start:mh1.end,'-','color','b','linewidth',2);
    % %             pause
    %             allval=ppval(mh2,mh2.start:mh2.end);
    %             allx=allval(1,:);ally=allval(2,:);
    %             plot3(allx,ally,mh2.start:mh2.end,'-','color','k','linewidth',2);
    % %             pause
    
    
    %             allval=ppval(sfit,t(1):t(end));
    %             allx=allval(1,:);ally=allval(2,:);
    %             plot3(allx,ally,t(1):t(end),'--','color','b','linewidth',3);
    %
    
    mh1=mhs(merspl(nm,2)); mh2=mhs(merspl(nm,1));
    tt1=mh1.start:mintime-xi; tt2=mintime+xi:mh2.end;
    
    % fit through splines
    xt1=ppval(mh1,tt1); xt2=ppval(mh2,tt2);
    xy=[xt1 xt2];
    t=[tt1 tt2];
    tr=t+0.01*rand(1,length(t))-0.005;
    addedMerged=addedMerged+1;
    sfit=splinefit(tr,xy,round(max(opt.minCPs,(tr(end)-tr(1))*opt.ncpsPerFrame )),fitnorm);
    sfit=adjustSplineStruct(sfit, min(t), max(t), alldpoints, T, 0, [], [], []);
    mhsnew(addedMerged) = sfit;
    
    
end
