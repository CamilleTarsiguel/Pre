function TNeighbors=getTemporalNeighbors(allpoints)
% get spatio-temporal neighbors for each detection 
% Es (Fig. 4 (a), PAMI)




global opt
nPoints=length(allpoints.xp);
TNeighbors=sparse(nPoints,nPoints);

allt=unique(sort(allpoints.tp));
allt=allt(1:end-1);
for t=allt
    thist=find(allpoints.tp==t);
    nextt=find(allpoints.tp==t+1);
    nt=numel(thist);
    ntt=numel(nextt);
%     if nt>1
        for np1=1:nt
            for np2=1:ntt
                eucldist=norm([allpoints.xp(thist(np1)) allpoints.yp(thist(np1))] - ...
                    [allpoints.xp(nextt(np2)) allpoints.yp(nextt(np2))]);
                if eucldist < opt.tau
                    TNeighbors(thist(np1),nextt(np2))=1;
                end                
            end
            
        end
%     end
    
end

% allt=unique(sort(allpoints.tp));
% allt=allt(1:end-2);
% for t=allt
%     thist=find(allpoints.tp==t);
%     nextt=find(allpoints.tp==t+2);
%     nt=numel(thist);
%     ntt=numel(nextt);
% %     if nt>1
%         for np1=1:nt
%             for np2=1:ntt
%                 eucldist=norm([allpoints.xp(thist(np1)) allpoints.yp(thist(np1))] - ...
%                     [allpoints.xp(nextt(np2)) allpoints.yp(nextt(np2))]);
%                 
%                 if eucldist < opt.tau*2
%                     TNeighbors(thist(np1),nextt(np2))=1;
%                 end                
%             end
%             
%         end
% %     end
%     
% end
% 
% allt=unique(sort(allpoints.tp));
% allt=allt(1:end-3);
% for t=allt
%     thist=find(allpoints.tp==t);
%     nextt=find(allpoints.tp==t+3);
%     nt=numel(thist);
%     ntt=numel(nextt);
% %     if nt>1
%         for np1=1:nt
%             for np2=1:ntt
%                 eucldist=norm([allpoints.xp(thist(np1)) allpoints.yp(thist(np1))] - ...
%                     [allpoints.xp(nextt(np2)) allpoints.yp(nextt(np2))]);
%                 
%                 if eucldist < opt.tau*3
%                     TNeighbors(thist(np1),nextt(np2))=1;
%                 end                
%             end
%             
%         end
% %     end
%     
% end

% TNeighbors=1*TNeighbors;
% TNeighbors=TNeighbors+TNeighbors';

end