function drawSplines(mh,used,labeling,allpoints,frames,lst,lw)
% plot splines color coded or gray if outliers
%


if nargin<6, lst='-'; end
if nargin<7, lw=2; end

global mhsglob
mhsglob = mh;

% show head and tail
head=0;
tail=0;
linesmoothing='on';
linesmoothing='off';
FontSize=12;
for id=used
    drawwidth=lw;
    mhs=mh(id);
    if labeling==0
        tt=linspace(1-tail,frames(end)+head,length(frames)*5);
        tt=linspace(mhs.start-tail,mhs.end+head,(mhs.end-mhs.start)*5);
        col=.6*ones(1,3);
        %         drawwidth=100*lw/splineGoodness(id);
        %         drawwidth=max(drawwidth,.1);
        %         drawwidth=min(drawwidth,3);
        %         col=500/splineGoodness(id)*ones(1,3); col(col>1)=1; col(col<0)=0; col=1-col;
        %         splineGoodness(id)
        %         drawwidth
        %         pause
    else
        supportPtsTs=allpoints.tp(labeling==id);
        supportPtsTs=sort(supportPtsTs);
        %         labeling
        %         id
        %         supportPtsTs
        %         [mh(id).start mh(id).end]
        tt=linspace(supportPtsTs(1)-tail,supportPtsTs(end)+head,length(supportPtsTs)*5);
        col=getColorFromID(id);
    end
    tt=mh(id).start:.5:mh(id).end;
    allval=ppval(mh(id),tt);
    allx=allval(1,:);ally=allval(2,:);
    
    %             col=getColorFromID(id);
    
    %     if tt(end)-tt(1)>30
    % find(labeling==id)
    % tt
    
    plot3(allx,ally,tt,lst,'color',col,'linewidth',drawwidth,'LineSmoothing',linesmoothing);
    text(allx(1),ally(1),tt(1),sprintf('%i',id),'color',col,'FontSize',FontSize);
    text(allx(end),ally(end),tt(end),sprintf('%i',id),'color',col,'FontSize',FontSize);
    %     end
    
    
    %% slope and curvature
    %     diff1=ppdiff(mh(id),1);
    %     diff2=ppdiff(mh(id),2);
    %     curvature=ppval(diff2,tt);
    %     slope=ppval(diff1,tt);
    %     abscurv=abs(curvature);sumabscurv=sum(abscurv);
    %     absslope=abs(slope);sumabsslope=sum(absslope);
    %     speed=sqrt(sum(slope.^2));
    %     sumabscurv=sumabscurv/max(sumabscurv);
    %     [slowest slowestt]=min(sumabsslope);
    %     normalizedsumabsslope=sumabsslope/max(sumabsslope);
    %     for np=1:length(tt)
    %         plot3(allx(np),ally(np),tt(np),'o','color',col,'MarkerSize',sumabsslope(np)/20);
    %     end
    %     plot3(allx(slowestt),ally(slowestt),tt(slowestt),'x','color',col,'MarkerSize',20);
    %     text(allx(slowestt),ally(slowestt),tt(slowestt),sprintf('%.1f',speed(slowestt)),'color',col);
    
    %     id\
    %     sl=getSplineSlope(mhs)
    %     [goodness components]=getSplineGoodness(mhs,1,allpoints,length(frames));
    %     fprintf('%10s%10s%10s%10s%10s%10s\n', ...
    %         'labelcost','fidlt','p_start','p_end','curvature','slope');
    %     fprintf('%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n', ...
    %         goodness,components(1),components(5),components(6),components(7),components(11));
    %
    %     pause
end
%     pause

head=0;
tail=0;

if tail
    % tail
    for id=used
        if labeling==0
            tt=linspace(1-tail,1,tail*10);
            col=.6*ones(1,3);
        else
            supportPtsTs=allpoints.tp(find(labeling==id));
            supportPtsTs=sort(supportPtsTs);
            tt=linspace(supportPtsTs(1)-tail,supportPtsTs(1),tail*10);
            tt=linspace(mh(id).start-tail,mh(id).start,tail*10);
            col=getColorFromID(id);
        end
        allval=ppval(mh(id),tt);
        allx=allval(1,:);ally=allval(2,:);
        plot3(allx,ally,tt,'-','color',min(col+.6,.8),'linewidth',1,'LineSmoothing','off');
        %     plot3(allx,ally,tt,'-','color',col,'linewidth',1,'LineSmoothing',linesmoothing);
        %     text(allx(1),ally(1),tt(1),sprintf('%i',id),'color',col);
    end
    
end

if head
    % head
    for id=used
        if labeling==0
            tt=linspace(frames(end),frames(end)+head,head*10);
            col=.6*ones(1,3);
        else
            supportPtsTs=allpoints.tp(find(labeling==id));
            supportPtsTs=sort(supportPtsTs);
            tt=linspace(supportPtsTs(end),supportPtsTs(end)+head,head*10);
            tt=linspace(mh(id).end,mh(id).end+head,head*10);
            col=getColorFromID(id);
        end
        allval=ppval(mh(id),tt);
        allx=allval(1,:);ally=allval(2,:);
        plot3(allx,ally,tt,'-','color',min(col+.7,.9),'linewidth',1,'LineSmoothing','off');
        %     plot3(allx,ally,tt,'-','color',col,'linewidth',1,'LineSmoothing',linesmoothing);
        %     text(allx(1),ally(1),tt(1),sprintf('%i',id),'color',col);
    end
    
end

drawnow