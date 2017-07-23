function [goodness, components]=getSplineGoodness(mh,used,alldpoints,T)
% compute the goodness (label cost) of each trajectory
% here based on spline structs with no derivatives

goodness=zeros(1,length(used));
components=zeros(length(used),5);
global sceneInfo opt

global occmaps

gF=opt.goodnessFactor;

% no need to compute if coefficient is 0 anyway
if ~gF, return; end


tau=sceneInfo.targetSize/2;


tt=1:T;
k=opt.conOpt.enParEfid(3);
ksq=k*k;

% convert splines to matrix representation
stateInfo=getStateFromSplines(mh, struct('F',T));

% go through all active ones (T*)
for m=used
    confpermodel=[];
	
    spl=mh(m);					% get current spline
    tt=spl.start:spl.end;		% evaluate from s to e
    splinexy=[stateInfo.X(:,m)';stateInfo.Y(:,m)'];
    
    
    detsWithin=-1*ones(1,T); % detections close to model
    detsWithinSig=-1*ones(1,T); % detections close to model
    oneDetWithinSig=-1*ones(1,T); % detections close to model

    insideAL=zeros(1,T); % inside areaLimits

    
    %% data fidelity (deprecated, dropped for PAMI paper)
    for t=tt
        sx=splinexy(1,t); sy=splinexy(2,t);
        
        
        detind=find(alldpoints.tp==t);
        insideAL(t)=1;
        detsWithin(t)=0;
        detsWithinSig(t)=0;
        oneDetWithinSig(t)=0;
   
        sp=[sx*ones(1,numel(detind)); sy*ones(1,numel(detind))];
        dp=[alldpoints.xp(detind);alldpoints.yp(detind)];
        d=sp-dp;
        
        alldist=sqrt(sum(d.^2));
        alldisthu=sqrt(sum(d.^2)./ksq + 1);

        
        
        inthisframe= alldist < tau;
%         inthisframe
        mindist=min(alldist);
        if isempty(detind), mindist=Inf; end
%         mindist
        siga=1;
        sigb=(tau)*siga;
        notinthisframeSig = 1-1./(1+exp(-siga*alldist+sigb));
        inthisframeSig = 1./(1+exp(-siga*alldist+sigb));
        inthisframeSig = 1./(1+exp(-siga*k*(alldisthu-1)+sigb));
        
%         alldist
%         inthisframe
%         prod(inthisframeSig)
%         pause
        
        
        confinframe=alldpoints.sp(detind(inthisframe>0));
        confpermodel=[confpermodel confinframe];
        
        if opt.occ
%         t
%         inthisframe

        
%         occmap=getDetOcclusionShadows(t);
%         occmap=occmaps(:,:,t);
%           eval(sprintf('occmap=occmap%i;',t));
        occmap=occmaps(t).map;
          
%         occmap2=occmap;
%         occmap2(round(syi),round(sxi))=0;
%         imshow(occmap2)
%         pause
        posx=max(1,round(syi(t))); posx=min(posx,sceneInfo.imgHeight);
        posy=max(1,round(sxi(t))); posy=min(posy,sceneInfo.imgWidth);        
        if occmap(posx,posy)
            inthisframe=[inthisframe 1];
        end
%         inthisframe
%         pause
        end
        
        detsWithin(t)=sum(inthisframe);
        detsWithinSig(t)=sum(inthisframeSig);
        oneDetWithinSig(t)=prod(inthisframeSig);
    end
    
    
    
    %% occlusion gaps
    nodetframes=[0 ~detsWithin(insideAL>0) 0];

   
    k=1;
    tbw=detsWithinSig(insideAL>0);
    allxgrkind= tbw>k; % always >=0

    tbw=(k^2/6)*(1-(1-((tbw)./k).^2).^3);
    tbw(allxgrkind)=k^2/6;
    tbw=6/k^2*tbw;
    
%     tbw
%     ~detsWithin(insideAL>0)
%     1-tbw
%     sum(abs(~detsWithin(insideAL>0)-(1-tbw)))
%     pause
    
    c=find(diff(nodetframes)==1); d=find(diff(nodetframes)==-1);
    occgaps=(d-c);
%     nodetframes
%     occgaps
%     [c;d]
    
    
    occgapspen=0;
    apppen=0;

    occgapsfac=opt.fidelityFactor;
%     occgaps
%     -log(1./(occgaps.^2))
%     pause
    switch(opt.fidFunction)
        case 1
            % original formulation (CVPR 2012)
            occgapspen=occgapsfac*sum(occgaps.^3);% + ...
        case 2
            % linear penalty
            occgapspen=occgapsfac*sum(occgaps);% + ...            
        case 3
            % lorentzian
            occgapspen=occgapsfac*sum(-log(1./(occgaps.^2)));
        case 4
            % new fid (CVPR 2013)
            occgapspen=occgapsfac*sum(oneDetWithinSig(insideAL>0));
    end
            

    notonedets=0;

    outsideAL=0;
    
    %% disttodets
    distfac=10;
    disttodets=0;
    
    spl=mh(m);
    
    
    %% persistence
    perfac=opt.persistenceFactor;
    persistencestart=0;
    persistenceend=0;
    
    if perfac
    

    s=spl.start; e=spl.end;
    %     xy=ppval(mh(m),[s e]);
    xy=[stateInfo.X([s e],m)';stateInfo.Y([s e],m)'];
    
    % from starting point
    if s>1.5
        xs=xy(1,1); ys=xy(2,1);
        if opt.track3d
            [dms,distances,~,~]=min_dist_im(xs,ys,sceneInfo.imOnGP);
        else
            [dms,distances,~,~]=min_dist_im(xs,ys,sceneInfo.imOnGP);
%             [dms,~,~]=min_distances(xs,ys,sceneInfo.trackingArea);
        end
        
        siga=5;
        sigb=(opt.borderMargin)*siga;
        sigscale=perfac*opt.borderMargin;
        
%         dms_=min(abs(dms),opt.borderMargin);
%         persistencestart=perfac*dms_
%         dms__=min(abs(dms/1000), opt.borderMargin/1000);
%         persistencestart=1000*perfac*dms__
%         pause
        
        

        k=opt.borderMargin/1000;
        y=(k^2/6)*(1-(1-((distances)./k).^2).^3);
        distgrthr=abs(distances)>k;
        y(distgrthr)=k^2/6;
        
        y=6/k^2*y;        y=y*k;
        persistencestart=prod(y)*perfac;
%         pause

    end
    
    % from end point
    if e<T-0.5
        xe=xy(1,2); ye=xy(2,2);
        if opt.track3d
            [dme,distances,~,~]=min_dist_im(xe,ye,sceneInfo.imOnGP);
        else
%             [dme,~,~]=min_distances(xe,ye,sceneInfo.trackingArea);
            [dme,distances,~,~]=min_dist_im(xe,ye,sceneInfo.imOnGP);
        end
% % %         dme=min(dme,opt.borderMargin);
% % %         persistenceend=perfac*dme;
% % %         
% % %         k=opt.borderMargin/1000;
% % %         y=k*(1-(1-((distances)./k).^2).^3);
% % %         distgrthr=abs(distances)>k;
% % %         y(distgrthr)=k;
% % % 
% % %         persistenceend=prod(y)*perfac;
% % % 
% % % %         persistenceend=0;
        siga=5;
        sigb=(opt.borderMargin)*siga;
        sigscale=perfac*opt.borderMargin;
        


        k=opt.borderMargin/1000;
        y=(k^2/6)*(1-(1-((distances)./k).^2).^3);
        distgrthr=abs(distances)>k;
        y(distgrthr)=k^2/6;
        
        y=6/k^2*y;        y=y*k;
        persistenceend=prod(y)*perfac;
        

    end

    end
    persistence = persistencestart + persistenceend;    
    
    
    %% dynamics 
    speedfac=1.;
    hang=0;

    nf=1; if opt.track3d, nf=1000; end

    switch(opt.speedFunction)
    
        case 1
	    %% dynamics curvature (original, CVPR 2012)
            hang=speedfac*max(abs(spl.coefs(:,1)));    
        case 2
            %% dynamics curvature (new)
            kappa=getSplineCurvature(spl,tt)*sceneInfo.frameRate;
    
            curvpenfac=opt.curvatureFactor;
            hang=sum(kappa)*curvpenfac;
        case 3
            %% dynamics curvature (new, CVPR 2013, PAMI) lorentzian
            kappa=getSplineCurvature(spl,tt)*sceneInfo.frameRate;
    
            curvpenfac=opt.curvatureFactor;
%             kappa
            hang=sum(-log(1./(kappa.^2+1)))*curvpenfac;
%             hang
%             pause
            
            if ~opt.track3d
                ypos=ppval(spl,linspace(tt(1),tt(end),length(tt)));
                ypos=ypos(2,:);
                scf=1;
                if sceneInfo.scenario==51 || sceneInfo.scenario==53
                    scf=min(1e5, 1./abs(ypos-200));
                end

                kappa=kappa .* scf; % dirty hack to scale dynamics in 2D
                hang=sum(-log(1./(kappa.^2+1)))*curvpenfac;
            end
    end
    
   
    
    %% penalize high speeds with slope sigmoid
    hlin=0;
    
    if opt.slopeFactor
    sigA=1/opt.tau;sigB=opt.tau*sigA;  
    sigA=0.02; sigB=300*sigA;
    sl=getSplineSlope(spl,tt)/nf*sceneInfo.frameRate;
    sp=opt.slopeFactor;

%     sl
    hlin=sp*sum(((sl-1).^2));
    if ~opt.track3d
        ypos=splinexy(2,tt);
%         
%         sl
%         ypos
%         ( 1./(ypos-200))
%         sl .* ( 1./(ypos-200))

        scf=1;
        if sceneInfo.scenario==51 || sceneInfo.scenario==53
            scf=min(1e5, 1./abs(ypos-200));% hack to scale dynamics in 2D (ETH)         
        end

        
        hlin=sp*sum( (sl .* scf).^2);
    end
    end
%     hlin=hlin+sp*min(((sl-1000).^2));

    % seg consistency
    [~, spl.index] = histc(tt,[-inf,spl.breaks(2:spl.pieces),inf]);
    stateVec=splinesToVec(spl);
    segpen = opt.segFactor*Eseg(stateVec,spl,opt);
    

    allpens=gF*[hlin hang persistence occgapspen segpen];
    components(m,:)=allpens;
    %  allpens
    goodness(m)=sum(allpens);
end
%% d(p,L) = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) )

goodness=goodness(used);
goodness(goodness>1e5)=1e5; % Integer overflow
end

function [dm, distances, ci, si]=min_dist_im(x,y,imOnGP)
% left
x0=imOnGP(1);y0=imOnGP(2);
x1=imOnGP(3);y1=imOnGP(4);
dl = p2l(x,y,x0,x1,y0,y1);

% top
x0=imOnGP(3);y0=imOnGP(4);
x1=imOnGP(5);y1=imOnGP(6);
du = p2l(x,y,x0,x1,y0,y1);

% right
x0=imOnGP(5);y0=imOnGP(6);
x1=imOnGP(7);y1=imOnGP(8);
dr = p2l(x,y,x0,x1,y0,y1);

% bottom
x0=imOnGP(7);y0=imOnGP(8);
x1=imOnGP(1);y1=imOnGP(2);
dd = p2l(x,y,x0,x1,y0,y1);
distances=[dl dr du dd];
% distances=abs(distances); % absolute
% distances=sqrt(1+(distances.^2))-1; % psudo huber




% choose the closest one
[dm, ci]=min(distances);
distances=distances/1000;

si=1;
% if x<minX || x>maxX || y<minY || y>maxY
%     si=-1;
% end
end

function dl=p2l(x,y,x0,x1,y0,y1)
% TODO, check for correctness
dl = ((y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) );
end

function [dm,  ci, si]=min_distances(x,y,areaLimits)
% returns min distance from x,y to border and index ci
error('take care of Tukey');
minX=areaLimits(1); % left border
maxX=areaLimits(2); % right border
minY=areaLimits(3); % bottom border
maxY=areaLimits(4); % top border

% determine distance to all four borders

% dist left
dl=abs(minX-x); % dl=x-minX;
% dist right
dr=abs(maxX-x); %dr=maxX-x;
% dist up
du=abs(minY-y); %du=y-minY;
% dist down
dd=abs(maxY-y); %dd=maxY-y;
distances=[dl dr du dd];

% choose the closest one
[dm, ci]=min(distances);

si=1;
if x<minX || x>maxX || y<minY || y>maxY
    si=-1;
end

end