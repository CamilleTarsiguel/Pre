function SNeighbors=getSpatialNeighbors(allpoints)
% Get spatial neighbors for each frame
% Ex in eq. (9), PAMI



global opt;
nPoints=length(allpoints.xp);
SNeighbors=sparse(nPoints,nPoints);

pts=[allpoints.xp; allpoints.yp];
allt=unique(sort(allpoints.tp));
for t=allt
    thist=find(allpoints.tp==t);
%         t
%         thist
    
    if length(thist)>1
        
        allpairs=nchoosek(thist,2);
%         allpairs

        % do not need extremely far ones
        nrm=zeros(1,size(allpairs,1));
        for ii=1:size(allpairs,1)
            nrm(ii)=norm(pts(:,allpairs(ii,1))-pts(:,allpairs(ii,2)));
        end
%         nrm
%         allpairs
        

		% size(allpairs)
		% nrm
		% nrm<opt.tau*10
		% nrm>opt.tau/2
        % TODO, move this pruning to pruneGraph.m
% 		tokeep=find(nrm<opt.tau*10 & nrm>opt.tau/2);
        tokeep=find(nrm>opt.tau/2);
        allpairs=allpairs(tokeep,:);
		% allpairs=allpairs(nrm>opt.tau/2,:);
%                 allpairs
%                 size(allpairs,1)
%                 pause

        if size(allpairs,1)
%             choosenew=nchoosek(1:size(allpairs,1),2);
%             choosenew
%             allpairs=allpairs(choosenew,:);
%             allpairs
%             sub2ind([nPoints nPoints],allpairs(:,1),allpairs(:,2))
            SNeighbors(sub2ind([nPoints nPoints],allpairs(:,1),allpairs(:,2)))=1;
%             pause
        end
    end
%     pause
%     thist
%     allpairs
%     pause
%     allpairs(:,1)
%     allpairs(:,2)
%     SNeighbors(allpairs(:,1),allpairs(:,2))=1;
%     (SNeighbors(1:10,1:10))
%     nextt=find(allpoints.tp==t+1);
%     nt=numel(thist);
%     ntt=numel(nextt);
%     if nt>1
%         for np1=1:nt
%             for np2=1:ntt
%                 eucldist=norm([allpoints.xp(thist(np1)) allpoints.yp(thist(np1))] - ...
%                     [allpoints.xp(nextt(np2)) allpoints.yp(nextt(np2))]);
%                 
%                 if eucldist < opt.tau
%                     TNeighbors(thist(np1),nextt(np2))=1;
%                 end                
%             end
%             
%         end
%     end
    
end

% allt=unique(sort(allpoints.tp));
% allt=allt(1:end-2);
% for t=allt
%     thist=find(allpoints.tp==t);
%     nextt=find(allpoints.tp==t+2);
%     nt=numel(thist);
%     ntt=numel(nextt);
%     if nt>1
%         for np1=1:nt
%             for np2=1:ntt
%                 eucldist=norm([allpoints.xp(thist(np1)) allpoints.yp(thist(np1))] - ...
%                     [allpoints.xp(nextt(np2)) allpoints.yp(nextt(np2))])/1000;
%                 
%                 if eucldist < .7
%                     TNeighbors(thist(np1),nextt(np2))=1;
%                 end                
%             end
%             
%         end
%     end
%     
% end

% TNeighbors=1*TNeighbors;
% TNeighbors=TNeighbors+TNeighbors';

end