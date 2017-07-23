function drawPoints(allpoints,labeling,outlierLabel,TNeighbors,Triplets)
% plot detections as points
% 


% return;

        global allflows
        global SNeighbors

npts=length(allpoints.xp);
for np=1:min(2000,npts) 
    msize=15;
    if length(labeling)>1,
        IDcol=getColorFromID(labeling(np)); marker='.';
        if labeling(np)==outlierLabel, IDcol='k'; marker='o'; msize=msize/2; end
    else
        IDcol=.6*ones(1,3); marker='.';
    end
%     plot3(allpoints.xp(np),allpoints.yp(np), allpoints.tp(np),marker,'color',IDcol,'MarkerSize',5);
            plot3(allpoints.xp(np),allpoints.yp(np), allpoints.tp(np),marker,'color',IDcol,'MarkerSize',msize*allpoints.sp(np));
    
            % show temporal pairs
%             thisneighb=find(TNeighbors(np,:));
%             nthisn=numel(thisneighb);
%             for nn=1:nthisn
%                 plot3([allpoints.xp(np) allpoints.xp(thisneighb(nn))], ...
%                     [allpoints.yp(np) allpoints.yp(thisneighb(nn))], ...
%                     [allpoints.tp(np) allpoints.tp(thisneighb(nn))],'color',[.5 .3 .8])
%             end
%             pause

        % show temporal triplets
%             triplets=find(Triplets(:,1)==np)';
%             nthisn=numel(triplets);
%             for nn=1:nthisn
%                 plot3(allpoints.xp(Triplets(triplets(nn),:)), ...
%                     allpoints.yp(Triplets(triplets(nn),:)), ...
%                 allpoints.tp(Triplets(triplets(nn),:)), ...
%                 'color',.6*ones(1,3))
%             end

% show spatial neighbors
%             thisneighb=find(SNeighbors(np,:));
%             nthisn=numel(thisneighb);
%             for nn=1:nthisn
%                 plot3([allpoints.xp(np) allpoints.xp(thisneighb(nn))], ...
%                     [allpoints.yp(np) allpoints.yp(thisneighb(nn))], ...
%                     [allpoints.tp(np) allpoints.tp(thisneighb(nn))],'color',[.5 .3 .8])
%             end
%         text(allpoints.xp(np),allpoints.yp(np), allpoints.tp(np),sprintf('%d',np));
%         pause

% show flow

%         plot3([allpoints.xp(np) allpoints.xp(np)+allflows.vx(np)], ...
%                     [allpoints.yp(np) allpoints.yp(np)+allflows.vy(np)], ...
%                     [allpoints.tp(np) allpoints.tp(np)+1],'color',[.5 .3 .8])
%             
% %         text(allpoints.xp(np),allpoints.yp(np), allpoints.tp(np),sprintf('%d',np));
%         pause

% show ori
%         plot3([allpoints.xp(np) allpoints.xp(np)+500*allpoints.dirx(np)], ...
%                     [allpoints.yp(np) allpoints.yp(np)+500*allpoints.diry(np)], ...
%                     [allpoints.tp(np) allpoints.tp(np)+1],'color',[.5 .3 .8])
%             
%         text(allpoints.xp(np),allpoints.yp(np), allpoints.tp(np),sprintf('%d',np));
%         pause
end

%% display pt info
randpts=randperm(npts); randpts=randpts(1:min(300,npts));
for np=randpts
    % pt number
%     text(allpoints.xp(np)-10,allpoints.yp(np), allpoints.tp(np),sprintf('%d',np));

%     % pt id
%     labid=0;
%     if length(labeling)>1, labid=labeling(np);
%         IDcol=getColorFromID(labeling(np)); 
%         if labeling(np)==outlierLabel, IDcol='k';  msize=msize/4; end
%     else
%         IDcol=.6*ones(1,3); 
%     end     
%     text(allpoints.xp(np),allpoints.yp(np), allpoints.tp(np),sprintf('%d',labid),'color',IDcol);
%     
    % pt score
%     text(allpoints.xp(np),allpoints.yp(np), allpoints.tp(np),sprintf('%.2f',allpoints.sp(np)));
end

drawnow;
end
