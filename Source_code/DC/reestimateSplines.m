function splines= ...
    reestimateSplines(allpoints,Nhood,used,labeling,nLabels,mhsall,DcostAll,opt,sceneInfo,flagOGM)
% refit splines to data points
% minimize continuous variables T
% 

if isempty(used)
    splines=mhsall;
    return
end

% TNeighbors=Nhood.TNeighbors;
% SNeighbors=Nhood.SNeighbors;

% if nargin<8, flagOGM=0; end


% global opt
T=opt.seqLength;

fitnorm='s';
if flagOGM ==1, fitnorm='a'; end
% fitnorm='s';
flagOGM=1;

nPoints=length(allpoints.xp);
minCPs=opt.minCPs;
ncpsPerFrame=opt.ncpsPerFrame;

exttail=2;
exthead=2;
nUsed=length(used);
Eold=zeros(0,4);Enew=zeros(0,4);
splines=mhsall;

F=length(used);
XGlob=zeros(T,F); 
YGlob=zeros(T,F); 
mcnt=0;
for m1=used
    mcnt=mcnt+1;
    tt=splines(m1).start:splines(m1).end;
    splinexy=ppval(splines(m1),tt);
    XGlob(tt,mcnt)=splinexy(1,:)'; 
    YGlob(tt,mcnt)=splinexy(2,:)'; 
end
stInfo.X=XGlob; stInfo.Y=YGlob;

% assert(flagOGM==1,'discrete optimization other than OpenGM is disabled');

if flagOGM
b=labelcount(labeling,nLabels);
labelUsed=~~b;
lu= labelUsed(1:end-1);
[proxGlob, proxtGlob, proxcostGlob]=getSplineProximity(splines(lu),T,stInfo,opt);
end
        
mcnt=0;
for m=used
    mcnt=mcnt+1;
%     mhsallnew=mhsall;
    
    mhsold=splines(m);
    
    %%%  Consistency check!
    % compute old energy
    Dcost=DcostAll;
    Lcost = getLabelCost(splines,opt);    
    h=setupGCO(nPoints,nLabels,Dcost,Lcost,[],[]);
    GCO_SetLabeling(h,labeling);
    [E, D, S, L] = GCO_ComputeEnergy(h);
    GCO_Delete(h);
    if flagOGM
        energy=evaluateFG(labeling,Dcost,Lcost,[],0,[],0);  
%           energy2=evaluateEnergy( ...
%               allpoints, Nhood,labeling, splines, opt, []);
%           energy
%           energy2
%          assert(isequal(energy,energy2));
        E=energy.value;
%         fprintf('old ens: '); fprintf('%d ',[Eogm Eun Epw Eex Ela Eprox]); fprintf('\n');
        

%         [prox proxt proxcost]=getSplineProximity(splines(lu),T, stInfo);
%         prox        
%         proxt
%         pause        
%         mcnt
        [prox2, proxt2, proxcost2]=getOneSplineProximity(splines(lu),T, stInfo, mcnt,opt,sceneInfo);
%         proxcost2
%         pause
        
%         prox2
%         proxt2
        proxcostNew=proxcostGlob;
        if mcnt==1
            proxcostNew(1,2:end)=proxcost2(2:end);
        elseif mcnt==length(used)
            proxcostNew(1:end-1,mcnt)=proxcost2(1:end-1)';
        else
            proxcostNew(1:mcnt-1,mcnt)=proxcost2(1:mcnt-1);
            proxcostNew(mcnt,mcnt+1:end)=proxcost2(mcnt+1:end);
        end

%         mcnt
%         [u v]=find(proxcost);
%         [u';v']
%         proxcost(~~proxcost)'
%         
%         proxcost2
% 
%         [u v]=find(proxcostNew);
%         [u';v']
%         proxcostNew(~~proxcostNew)'
%         
%         sum(abs(proxcost(:)-proxcostNew(:)))
%         assert(sum(abs(proxcost(:)-proxcostNew(:)))<1e-5);
        proxcost=proxcostNew;
%         pause
    
    
%         [proxa proxta proxcosta]=getSplineProximity(splines(lu),T, stInfo);
%         tmpproxcosta=proxcosta+proxcosta';
%         tmpsuma=sum(tmpproxcosta(mcnt,:));

        tmpproxcost=proxcost+proxcost';
        tmpsum=sum(tmpproxcost(mcnt,:));
%         assert(abs(tmpsum-tmpsuma)<1e-5);
        
%         tmpsum=sum(proxcost(:));
%         fprintf('%d %d %d\n',E, tmpsum, E+tmpsum);
        E=E+tmpsum;
        Eogm=E;
        Eprox=tmpsum;
    end
    

%     Eold=[Eold; E D S L];
    
    % find all points labeled m
    supportPts=find(labeling==m);
    nSP=length(supportPts);
    
%     fprintf('from reestimatedSplines.m\n');
%     m
%     mhsold
%     labeling
%     supportPts
    % careful! if less than 2 points, make standing object
    if nSP<1
        error('This cannot happen... I think');
    elseif nSP==1
        xy=[allpoints.xp(supportPts);allpoints.yp(supportPts)];        xy=[xy xy];
        t=allpoints.tp(supportPts); 
        confs=allpoints.sp(supportPts); 
        if t==1,
            t=[t t+1];
        else
            t=[t-1 t];
        end
        torig=t;
        confsorig=confs;
        confs=[confs confs];

        supPtsFrames=t;
    else
        xy=[allpoints.xp(supportPts);allpoints.yp(supportPts)];
        t=allpoints.tp(supportPts);
        confs=allpoints.sp(supportPts); 
        [supPtsFrames, sortind]=sort(allpoints.tp(supportPts));    
               
        % try here this...
        % take n first and last points and fit a line to extrapolate
        % add these to support
        takenpts=4;
        torig=t;
        confsorig=confs;
        if nSP>=takenpts
            if exttail
                ntail=sortind(1:takenpts);     
                tailpts=supportPts(ntail);
                xytail=[allpoints.xp(tailpts);allpoints.yp(tailpts)];
                ttail=allpoints.tp(tailpts);
                ctail=allpoints.sp(tailpts);
            
                if length(unique(ttail))>2
                    tailline=splinefit(ttail,xytail,1,2,ctail,fitnorm);
                    taillinepts=ppval(tailline,supPtsFrames(1)-exttail:supPtsFrames(1)-1);
                    xy=[taillinepts xy];
                    t=[supPtsFrames(1)-exttail:supPtsFrames(1)-1 torig];
                    confs=[ones(1,exttail) confsorig];
                end
            end
            
            if exthead
                nhead=sortind(end-takenpts+1:end);
                headpts=supportPts(nhead);
                xyhead=[allpoints.xp(headpts);allpoints.yp(headpts)];
                thead=allpoints.tp(headpts);
                chead=allpoints.sp(headpts);
                if length(unique(thead))>2
                    headline=splinefit(thead,xyhead,1,2,chead,fitnorm);
                    headlinepts=ppval(headline,supPtsFrames(end)+1:supPtsFrames(end)+exthead);

                    xy=[xy headlinepts];
                    t=[t supPtsFrames(end)+1:supPtsFrames(end)+exthead];
                    confs=[confs ones(1,exthead)];
                end
            end
        end
    end
    
    splineorder=4;        
    order=min(nSP,splineorder);
    
    
    
    trackLength(m)=supPtsFrames(end)-supPtsFrames(1);            
    ncps=max(minCPs,round(trackLength(m)*ncpsPerFrame));
    tr=t;
%     tr=t+opt.randFit*rand(1,length(t)); % add random noise to avoid NaN in fitting (LOOK INTO THIS!)    
%     breaks=linspace(supPtsFrames(1),supPtsFrames(end),ncps);
    tryfit=splinefit(tr,xy,mhsold.pieces,mhsold.order,confs,fitnorm);


    % force cubic spline
    if order<splineorder
        sortedt=sort(tr);
        t=linspace(sortedt(1),sortedt(end),splineorder);
        xy=ppval(tryfit,t);
        sfit=splinefit(t,xy,1,splineorder,fitnorm);
    else
        sfit=tryfit;
    end
    sfit=adjustSplineStruct(sfit, min(torig), max(torig), allpoints, T, 0, [], [], []);

%         if m==43
%         xy
%         t
%         confs
%         sfit
%         sfit.coefs
%     end    
    % consistency check
    % compute new energy
    spl2=splines; spl2(m)=sfit;
    DcostNew=getUnarySpline(2,nPoints,sfit,allpoints,opt);
    Dcost=DcostAll;
    Dcost(m,:)=DcostNew(1,:);
    Lcost = getLabelCost(spl2,opt);
%     h=setupGCO(nPoints,nLabels,Dcost,Lcost,[],[]);
%     GCO_SetLabeling(h,labeling);
%     [E_, D_, S_, L_] = GCO_ComputeEnergy(h);    
    if flagOGM
        energy_=evaluateFG(labeling,Dcost,Lcost,[],0,[],0);    
        E_=energy_.value;
%         Eun_=Eun; Epw_=Epw; Eex_=Eex; Ela_=Ela; 
%         fprintf('new ens: '); fprintf('%d ',[Eogm Eun Epw Eex Ela Eprox]); fprintf('\n');

%         b=labelcount(labeling,nLabels);
%         labelUsed=~~b;
%         lu= labelUsed(1:end-1);
        stInfo2=stInfo;
        tt=sfit.start:sfit.end;
        splinexy=ppval(sfit,tt);
        stInfo2.X(tt,mcnt)=splinexy(1,:)'; 
        stInfo2.Y(tt,mcnt)=splinexy(2,:)'; 
        
%         [prox proxt proxcost]=getSplineProximity(spl2(lu),T, stInfo2);

        [prox2, proxt2, proxcost2]=getOneSplineProximity(spl2(lu),T, stInfo2,mcnt,opt,sceneInfo);
        
%         prox2
%         proxt2
        proxcostNew=proxcostGlob;
        if mcnt==1
            proxcostNew(1,2:end)=proxcost2(2:end);
        elseif mcnt==length(used)
            proxcostNew(1:end-1,mcnt)=proxcost2(1:end-1)';
        else
            proxcostNew(1:mcnt-1,mcnt)=proxcost2(1:mcnt-1);
            proxcostNew(mcnt,mcnt+1:end)=proxcost2(mcnt+1:end);
        end
        
%         proxcost
%         proxcost2
%         proxcostNew
%         proxcost-proxcostNew
%         assert(sum(abs(proxcost(:)-proxcostNew(:)))<1e-5);
        proxcost=proxcostNew;

        
%         [proxa proxta proxcosta]=getSplineProximity(spl2(lu),T, stInfo2);
%         tmpproxcost=proxcosta+proxcosta';
%         tmpsuma_=sum(tmpproxcost(mcnt,:));
%         proxcost
%         proxcosta
        
        tmpproxcost=proxcost+proxcost';
        tmpsum_=sum(tmpproxcost(mcnt,:));
%         tmpsum_
%         tmpsuma_
%         assert(abs(tmpsum_-tmpsuma_)<1e-5);
        
%         tmpsum_=sum(proxcost(:));
%         fprintf('%d %d %d\n',E_, tmpsum_, E_+tmpsum_);
        E_=E_+tmpsum_;
        Eprox_=tmpsum_;
        Eogm_=E_;
    end
    
%     GCO_Delete(h);
%     Enew=[Enew; E_ D_ S_ L_];
       
    % only if continuous optimization reduced
%     printLabelCost(mhsold)
%     printLabelCost(sfit)
%     figure(1)
%     drawSplines(sfit,1,0,allpoints,T);
    
    if E_<=E
%         energy
%         energy_
%         fprintf('Number %d refit, old: %f, new %f, diff: %f\n',m,E,E_,E-E_);
        splines(m)=sfit;
        DcostAll=Dcost;
		
		if flagOGM
            proxcostGlob=proxcost;
%         stInfo.X(:,mcnt)=0;stInfo.Y(:,mcnt)=0;
%         tt=sfit.start:sfit.end;
%         splinexy=ppval(sfit,tt);
%         stInfo.X(tt,mcnt)=splinexy(1,:)'; 
%         stInfo.Y(tt,mcnt)=splinexy(2,:)'; 
            stInfo=stInfo2;
		end

    else
%         fprintf('Number %d not refit, old: %f, new %f\n',m,E,E_)
    end
%     fprintf('%f\t',([Eogm Eun Epw Eex Ela Eprox])); fprintf('\n');
%     fprintf('%f\t',([Eogm_ Eun_ Epw_ Eex_ Ela_ Eprox_])); fprintf('\n');

   % if tmpcnt<20, save(sprintf('tmp/tmp%d/tmpvar_%05d.mat',processnr,tmpcnt),'*'); tmpcnt=tmpcnt+1; end

%     fprintf('\n');
    
%     pause
end

end