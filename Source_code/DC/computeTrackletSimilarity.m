%% create tracklet similarity matrix
stateInfo.targetsExist=getTracksLifeSpans(stateInfo.X);
N=size(stateInfo.X,2);
F=size(stateInfo.X,1);
speeds=Inf*ones(N);
gaps=Inf*ones(N);
distances=Inf*ones(N);
distances1=Inf*ones(N);
distances2=Inf*ones(N);
distS=Inf*ones(N);
distS1=Inf*ones(N);
distS2=Inf*ones(N);
distBhatt=Inf*ones(N);

% headsh=Inf



lw=3;
supptimelength=5;
simmatfile=sprintf('tmp/linked/simMat-clnr%03d-s%d.mat',clusternr,scenario);


if 1 && exist(simmatfile,'file')
    load(simmatfile);
else
    


% fit linear splines
% allheadlines=[];alltaillines=[];
% allheadlhines=[];alltailhlines=[];
clear allhead* alltail*
for m1=1:N
    if ~mod(m1,10), fprintf('.'); end

    tail1t=stateInfo.targetsExist(m1,1);
    head1t=stateInfo.targetsExist(m1,2);
    exfr1=tail1t:head1t;
    head1x=stateInfo.X(head1t,m1);
    head1y=stateInfo.Y(head1t,m1);
    head1h=stateInfo.H(head1t,m1);
    if length(exfr1)<supptimelength, continue; end

    % fit head
    headtimes=exfr1(end-supptimelength+1:end);
    headxy=[stateInfo.X(headtimes,m1)';stateInfo.Y(headtimes,m1)'];
    headline=splinefit(headtimes,headxy,1,2);
    allheadlines(m1)=headline;
    
    headh=stateInfo.H(headtimes,m1)';
    headhline=splinefit(headtimes,headh,1,2);
    allheadhlines(m1)=headhline;
    
    % fit tail
    tail1x=stateInfo.X(tail1t,m1);
    tail1y=stateInfo.Y(tail1t,m1);
    tail1h=stateInfo.H(tail1t,m1);
        
    tailtimes=exfr1(1:supptimelength);
    tailxy=[stateInfo.X(tailtimes,m1)';stateInfo.Y(tailtimes,m1)'];
    tailline=splinefit(tailtimes,tailxy,1,2);
    alltaillines(m1)=tailline;
        
    tailh=stateInfo.H(tailtimes,m1)';
    tailhline=splinefit(tailtimes,tailh,1,2);
    alltailhlines(m1)=tailhline;
end
fprintf('\n');
% pause
%%
fname=sprintf('tmp/allimages-s%d.mat',scenario);

evstr='size(allimages,4)==F';
if ~exist('allimages','var')
    if exist(fname,'file')
        load(fname);
    else
        allimages=zeros(sceneInfo.imgHeight,sceneInfo.imgWidth,3,length(sceneInfo.frameNums),'uint8');
        for t=1:length(sceneInfo.frameNums)
            if ~mod(t,10), fprintf('.'); end
            allimages(:,:,:,t)=imread([sceneInfo.imgFolder sprintf(sceneInfo.imgFileFormat,sceneInfo.frameNums(t))]);
        end
        save(fname,'allimages')
        fprintf('\n');
    end
elseif ~eval(evstr)
    load(fname);   
end

%%
% pause
for m1=1:N
    if ~mod(m1,10), fprintf('.'); end
           
    tail1t=stateInfo.targetsExist(m1,1);
    head1t=stateInfo.targetsExist(m1,2);
    exfr1=tail1t:head1t;
    head1x=stateInfo.X(head1t,m1);
    head1y=stateInfo.Y(head1t,m1);
    head1h=stateInfo.H(head1t,m1);
    
    if length(exfr1)<supptimelength, continue; end

    headline=allheadlines(m1);
    headhline=allheadhlines(m1);
    
    for m2=m1+1:N
        tail2t=stateInfo.targetsExist(m2,1);
        head2t=stateInfo.targetsExist(m2,2);
        exfr2=tail2t:head2t;
        
        if length(exfr2)<supptimelength, continue; end
        
        tail2x=stateInfo.X(tail2t,m2);
        tail2y=stateInfo.Y(tail2t,m2);
        tail2h=stateInfo.H(tail2t,m2);
        
        tailline=alltaillines(m2);
        tailhline=alltailhlines(m2);

        timegap=tail2t-head1t;
        gaps(m1,m2)=timegap;
        distances(m1,m2)=norm([head1x head1y]-[tail2x tail2y]);
        
        extrap1=ppval(headline,tail2t);
        extrap2=ppval(tailline,head1t);
        
        extrap1h=ppval(headhline,tail2t);
        extrap2h=ppval(tailhline,head1t);
        
        distances1(m1,m2)=norm(extrap1-[tail2x tail2y]')/tail2h;
        distances2(m1,m2)=norm(extrap2-[head1x head1y]')/head1h;
        distS(m1,m2)=abs(extrap1h-extrap2h);%/min(extrap1h,extrap2h);
        distS1(m1,m2)=abs(extrap1h-tail2h)/tail2h;
        distS2(m1,m2)=abs(extrap2h-head1h)/head1h;
        
        speed=distances(m1,m2)/timegap;
        speeds(m1,m2)=speed;
        
        if timegap>0 && timegap <50
            
%             head1im=double(imread([sceneInfo.imgFolder sprintf(sceneInfo.imgFileFormat,sceneInfo.frameNums(head1t))]))/255;
            head1im=double(allimages(:,:,:,head1t))/255;
            [head1th1 head1th2]=getTargetHist(head1im,stateInfo,head1t,m1);
%             pause
            
%             tail2im=double(imread([sceneInfo.imgFolder sprintf(sceneInfo.imgFileFormat,sceneInfo.frameNums(tail2t))]))/255;
            tail2im=double(allimages(:,:,:,tail2t))/255;

            [tail2th1 tail2th2]=getTargetHist(tail2im,stateInfo,tail2t,m2);


            mbb1=multibandBhattacharyya(head1th1',tail2th1');
            mbb2=multibandBhattacharyya(head1th2',tail2th2');
            distBhatt(m1,m2)=mbb1+mbb2;

%             title(sprintf('%.2f',mbb1))
%             pause

            

            
            
        end
    end
end
clear allimages
fprintf('\n');

sumdist=distances1+distances2;
sumdistS=distS1+distS2;
save(simmatfile,'sumdist','sumdistS','distances*','distS*','distBhatt','gaps');
end
% pause
%%

[metrics2d metrics3d addInfo2d addInfo3d]=printFinalEvaluation(stateInfo,gtInfo,sceneInfo,opt);
stInfo=stateInfo;
C=sumdist;
C=sumdistS;
% C=C+max(sumdist(sumdist<Inf))/max(distBhatt(distBhatt<Inf))*distBhatt + ...
%     max(sumdist(sumdist<Inf))/max(distS(distS<Inf))*distS;

maxDist=max(sumdist(sumdist<Inf).^2);


nrmlzDist=sumdist.^2;
nrmlzScale=sumdistS.^2;
nrmlzBhat=distBhatt;
nrmlzGaps=gaps;

% nrmlzDist=normToZeroOne(nrmlzDist);
% nrmlzScale=normToZeroOne(nrmlzScale);
% nrmlzBhat=normToZeroOne(nrmlzBhat);
% nrmlzGaps=normToZeroOne(nrmlzGaps);

featnormfile=sprintf('tmp/linked/simMat-clnr%03d-featnorm.mat',clusternr);
load(featnormfile);
nrmlzDist=(nrmlzDist-featmeanD)/stddevD;
nrmlzScale=(nrmlzScale-featmeanS)/stddevS;
nrmlzBhat=(nrmlzBhat-featmeanB)/stddevB;
nrmlzGaps=(nrmlzGaps-featmeanG)/stddevG;

wt.dist=1; % leave at one
wt.scal=1.2;
wt.appe=0.002;
wt.gaps=0.0001;


rng(exper);
if exper~=0
wt.scal=rand*2*wt.scal;
wt.appe=rand*2*wt.appe;
wt.gaps=rand*2*wt.gaps;
end

C= ...
    wt.dist * nrmlzDist+ ...
    wt.scal * nrmlzScale + ...
    wt.appe * nrmlzBhat + ...
    wt.gaps * nrmlzGaps;

% C=sumdistS;
numel(find(isnan(C)))
C(gaps<1)=Inf;
numel(find(isnan(C)))

%
firstmin=min(C(:));

firstMOTA=metrics2d(12);
remids=[];
allmind=[];
% pause

allm2d=zeros(0,14);
    [mind, bestpair]=min(C(:));
    allmind=[allmind mind];
    [b1, b2]=ind2sub(size(C),bestpair);
% while  9.170714*firstmin>min(C(:))
% while    1.843023*firstmin>min(C(:))
% while 3*firstmin>min(C(:))
% while min(C(:))<-0.429078
% while min(C(:))<-0.51562
% while min(C(:))<-0.70986
% while min(C(:))<-0.520493
% while min(C(:))< -0.520493
% while min(C(:))<  -0.398305
% while min(C(:))<-0.313602

% loopscale=C(b1,b2)/firstmin;
% while loopscale>0.998947 
% while min(C(:))<-0.313602
while min(C(:))<-0.452234
% while length(remids)<round(N/4) %&& firstMOTA<1.2*metrics2d(12)
    
    [mind, bestpair]=min(C(:));
    allmind=[allmind mind];
    [b1, b2]=ind2sub(size(C),bestpair);
%     [b1 b2]
    fprintf('sim: %f,  scale: %f\n', C(b1,b2), C(b1,b2)/firstmin);
%     1*sumdist(b1,b2)^2
    remids=[remids b2];
    m1=b1;m2=b2;
    length(remids)
    
%     exfr1=stateInfo.targetsExist(b1,1):stateInfo.targetsExist(b1,2);
%     exfr2=stateInfo.targetsExist(b2,1):stateInfo.targetsExist(b2,2);
    exfr1=find(stInfo.X(:,b1))';
    exfr2=find(stInfo.X(:,b2))';
    
%         figure(5)
%     clf
%     plot3(stInfo.X(exfr1,b1),stInfo.Y(exfr1,b1),exfr1, 'b', ...
%         stInfo.X(exfr2,b2),stInfo.Y(exfr2,b2),exfr2,'r','linewidth',lw)
%     xlim([1 640]);
%     ylim([240 480]);
%     box on
%     grid on
%     hold on
%     
%     
%     
%     headtimes=exfr1(end-supptimelength+1:end);
%     headxy=[stInfo.X(headtimes,m1)';stInfo.Y(headtimes,m1)'];
%     headline=splinefit(headtimes,headxy,1,2);
%     tailtimes=exfr2(1:supptimelength);
%     tailxy=[stInfo.X(tailtimes,m2)';stInfo.Y(tailtimes,m2)'];
%     tailline=splinefit(tailtimes,tailxy,1,2);
%     
%     
%     extraptimes=exfr1(end):exfr2(1);
%     headextrap=ppval(headline,extraptimes);
%     tailextrap=ppval(tailline,extraptimes);
%     plot3(headextrap(1,:),headextrap(2,:),extraptimes,'b', ...
%         tailextrap(1,:),tailextrap(2,:),extraptimes,'r')
    
    stInfoTMP=stInfo;
    stInfo=mergeTracks(stInfo, b1, b2);
    
%     plot3(stInfo.X(exfr1(1):exfr2(end),b1),stInfo.Y(exfr1(1):exfr2(end),b1),exfr1(1):exfr2(end), 'k')
    
    stInfEv=stInfo;
    keepids=setdiff(1:size(stInfo.X,2),remids);
    stInfEv.X=stInfEv.X(:,keepids);stInfEv.Y=stInfEv.Y(:,keepids);
    stInfEv.Xi=stInfEv.Xi(:,keepids);stInfEv.Yi=stInfEv.Yi(:,keepids);
    stInfEv.W=stInfEv.W(:,keepids);stInfEv.H=stInfEv.H(:,keepids);
%     [metrics2d metrics3d addInfo2d addInfo3d]=printFinalEvaluation(stInfEv,gtInfo,sceneInfo,opt);
    allm2d(end+1,:)=metrics2d;
    
%     pause
%     
%     figure(6)
%     headh=stInfo.H(headtimes,m1)';
%     headhline=splinefit(headtimes,headh,1,2);
% 
%     tailh=stInfo.H(tailtimes,m1)';
%     tailhline=splinefit(tailtimes,tailh,1,2);
%     
%     
%     headhextrap=ppval(headhline,extraptimes);
%     tailhextrap=ppval(tailhline,extraptimes);
%     
%     clf
%     hold on
%     plot(exfr1,stInfo.H(exfr1,b1), 'b', ...
%         exfr2,stInfo.H(exfr2,b1),'r','linewidth',lw)
% 
%     plot(extraptimes,headhextrap,'b',extraptimes,tailhextrap,'r');
% 
%     plot(exfr1(1):exfr2(end),stInfo.H(exfr1(1):exfr2(end),b1), 'k')
%     pause
%     figure(7)
%     clf
%     % track 1
%     for t=exfr1(end)%exfr1(end-supptimelength+1:end):exfr2(1:supptimelength)
%         
%         subplot(121); cla
%         im=double(imread([sceneInfo.imgFolder sprintf(sceneInfo.imgFileFormat,sceneInfo.frameNums(t))]))/255;
%         imshow(im,'Border','tight')
%         hold on
%         if stInfoTMP.Xi(t,b1)
%             bleft= stInfoTMP.Xi(t,b1)- stInfoTMP.W(t,b1)/2;            bright= stInfoTMP.Xi(t,b1)+ stInfoTMP.W(t,b1)/2;
%             btop= stInfoTMP.Yi(t,b1)- stInfoTMP.H(t,b1);            bbottom= stInfoTMP.Yi(t,b1);
%             boxcurvature=[.3,.3*( stInfoTMP.W(t,b1)/ stInfoTMP.H(t,b1))]; boxcurvature=max(0,boxcurvature);boxcurvature=min(1,boxcurvature);
%             rectangle('Position',[bleft,btop, stInfoTMP.W(t,b1), stInfoTMP.H(t,b1)],'Curvature',boxcurvature,'EdgeColor',getColorFromID(1),'linewidth',2);
%         end
%     end
%        
%     for t=exfr2(1)
%         subplot(122); cla
%         im=double(imread([sceneInfo.imgFolder sprintf(sceneInfo.imgFileFormat,sceneInfo.frameNums(t))]))/255;
%         imshow(im,'Border','tight')
%         hold on
%         if stInfoTMP.Xi(t,b2)
%             bleft= stInfoTMP.Xi(t,b2)- stInfoTMP.W(t,b2)/2;            bright= stInfoTMP.Xi(t,b2)+ stInfoTMP.W(t,b2)/2;
%             btop= stInfoTMP.Yi(t,b2)- stInfoTMP.H(t,b2);            bbottom= stInfoTMP.Yi(t,b2);
%             boxcurvature=[.3,.3*( stInfoTMP.W(t,b2)/ stInfoTMP.H(t,b2))]; boxcurvature=max(0,boxcurvature);boxcurvature=min(1,boxcurvature);
%             rectangle('Position',[bleft,btop, stInfoTMP.W(t,b2), stInfoTMP.H(t,b2)],'Curvature',boxcurvature,'EdgeColor',getColorFromID(1),'linewidth',2);
%         end
%         pause(.01)
%     end
%     pause

    

%     C(b1,b2)=Inf;    
    C(b1,:)=C(b2,:);
    C(b2,:)=Inf;
    C(:,b2)=Inf;

    [mind, bestpair]=min(C(:));
    allmind=[allmind mind];
    [b1, b2]=ind2sub(size(C),bestpair);
    loopscale=C(b1,b2)/firstmin;

%     C(b1,:)=Inf;
%     C(:,b2)=Inf;
%     C=C(keepids,:);
%     C=C(:,keepids);
    
end


%%
keepids=setdiff(1:N,remids);
stInfo.X=stInfo.X(:,keepids);stInfo.Y=stInfo.Y(:,keepids);
stInfo.Xi=stInfo.Xi(:,keepids);stInfo.Yi=stInfo.Yi(:,keepids);
% stInfo.Xgp=stInfo.Xgp(:,keepids);stInfo.Ygp=stInfo.Ygp(:,keepids);
stInfo.W=stInfo.W(:,keepids);stInfo.H=stInfo.H(:,keepids);

if savesim
    resfiledir=sprintf('tmp/linked/simres/%03d/',clusternr);
    if ~exist(resfiledir,'dir'), mkdir(resfiledir); end
    resfile=sprintf('%sexp%03d-s%02d.mat',resfiledir,exper,scenario);
    save(resfile,'stInfo','wt','allm2d','allmind','first*');
end