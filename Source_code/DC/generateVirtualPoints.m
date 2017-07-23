function [newPoints, newLabeling]=...
    generateVirtualPoints(allpoints,labeling, used)
% THIS IS NOT USED ANYMORE

% generate virtual points to fit the spline on the safety margin
% alldpoints - all old points
% labeling - current labeling, assuming outliers are -1

newPoints=allpoints;
newLabeling=labeling;
return;
% outlierLabel=-1;


% safety margin = 2 frames in either direction
exttail=2;
exthead=2;

for id=used
    supportPts=find(labeling==id);
    
    % How many unique time steps?
    uniqueT=unique(allpoints.tp(supportPts));
    
    % initialize with empty sets
    newxy=[]; newt=[];
    
    % determine locations and times
    if length(uniqueT) < 1
        warning('This cannot happen... I think');
    elseif length(uniqueT)==1
        % if only one point (or time step), extend at the same (mean) location
        t=allpoints.tp(supportPts(1));
        xy=[mean(allpoints.xp(supportPts));mean(allpoints.yp(supportPts))];
        newxy=repmat(xy,1,exttail+exthead);
        newt=[t-exttail:t-1 t+1:t+exthead];
    else
        % otherwise linear extrapolate to both ends
        
        % sort points according to their time stamps
        [supPtsFrames, sortind]=sort(allpoints.tp(supportPts));
        

                
        if length(uniqueT)<4
            % if two or three, fit through all those       
            xytail=[allpoints.xp(supportPts);allpoints.yp(supportPts)];
            ttail=allpoints.tp(supportPts);
            ctail=allpoints.sp(supportPts);
            
            % tail and head are the same in this case
            xyhead=xytail;            thead=ttail; chead=ctail;

        else
            % otherwise find all points with first four (resp. last four)
            % distinct frames and fit through those
            
            ntail= allpoints.tp(supportPts) <= uniqueT(4);
            tailpts=supportPts(ntail);
            xytail=[allpoints.xp(tailpts);allpoints.yp(tailpts)];
            ttail=allpoints.tp(tailpts);
            ctail=allpoints.sp(tailpts);
            
            nhead=allpoints.tp(supportPts) >= uniqueT(end-3);
            headpts=supportPts(nhead);
            xyhead=[allpoints.xp(headpts);allpoints.yp(headpts)];
            thead=allpoints.tp(headpts);
            chead=allpoints.sp(headpts);
%             [supPtsFrames, sortind]=sort(allpoints.tp(supportPts));
            
        end                
        
        % actual fit, and evaluate beyond start / end
        tailline=splinefit(ttail,xytail,1,2,ctail);
        newttail=supPtsFrames(1)-exttail:supPtsFrames(1)-1;        
        newxytail=ppval(tailline,newttail);
        
%         chead
        headline=splinefit(thead,xyhead,1,2, chead);
        newthead=supPtsFrames(end)+1:supPtsFrames(end)+exthead;
        newxyhead=ppval(headline,newthead);
       
        % finally, add new coordinates and times
        newxy=[newxytail newxyhead];
        newt=[newttail newthead];

    end
    
    
    newPoints.xp=[newPoints.xp newxy(1,:)];
    newPoints.yp=[newPoints.yp newxy(2,:)];
    newPoints.sp=[newPoints.sp ones(1,length(newt))];
    newPoints.tp=[newPoints.tp newt];
    newLabeling=[newLabeling id*ones(1,length(newt))];
    
end