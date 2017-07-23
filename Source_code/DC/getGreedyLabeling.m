function labeling=getGreedyLabeling(alldpoints, mh,T)


npts=length(alldpoints.xp);
stateInfo=getStateFromSplines(mh, struct('F',T));
stateInfo=postProcessState(stateInfo);
nLabels=length(mh)+1;

labeling=nLabels*ones(1,npts);
for np=1:npts
    t=alldpoints.tp(np);
    xd=alldpoints.xp(np);
    yd=alldpoints.yp(np);

    
    % find closest spline
    extar=find(stateInfo.X(t,:));
    
    if numel(extar)
    
        splpos=[stateInfo.X(t,extar);stateInfo.Y(t,extar)];
        nsplpos=size(splpos,2);
        
        reppt=repmat([xd;yd],1,nsplpos);
        ddists=sqrt(sum((splpos-reppt).^2));
            
            %         pause
            [mindist ddistsI]=min(ddists);
            if mindist<1000
                labeling(np)=ddistsI;
            end

    end
    
    
end

end