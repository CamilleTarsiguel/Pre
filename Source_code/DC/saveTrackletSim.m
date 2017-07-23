%% create and save tracklet similarity matrix (ETH only)
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

if ~exist(fname,'file')
allimages=zeros(sceneInfo.imgHeight,sceneInfo.imgWidth,3,length(sceneInfo.frameNums),'uint8');
for t=1:length(sceneInfo.frameNums)
    if ~mod(t,10), fprintf('.'); end
    allimages(:,:,:,t)=imread([sceneInfo.imgFolder sprintf(sceneInfo.imgFileFormat,sceneInfo.frameNums(t))]);
end
save(fname,'allimages')
fprintf('\n');
end
evstr='size(allimages,4)==F';
if ~exist('allimages','var')
    load(fname);
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