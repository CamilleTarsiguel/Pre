function mhs=getSplinesFromGT(Xgt,Ygt,frames,alldpoints,T)
% Fit splines through the ground truth trajectories
% FOR DEBUGGING ONLY

global opt

newmods=0;

mhs=getEmptyModelStruct;


[Fgt, Ngt]=size(Xgt);
    if length(frames)<Fgt
        Xgt=Xgt(frames,:); Ygt=Ygt(frames,:);
    end

exttail=5;
exthead=5;
for id=1:Ngt
    
    
    cptimes=find(Xgt(:,id));
    if length(cptimes)>3
        %% extend to first and last frame
        ff=cptimes(1);  lf=cptimes(end); trl=lf-ff;
        ncps=max(opt.minCPs,round(trl*opt.ncpsPerFrame));
        
        torig=cptimes';
        xy=[Xgt(cptimes,id) Ygt(cptimes,id)]';

        
        ttail=cptimes(1:4);
        xytail=[Xgt(ttail,id) Ygt(ttail,id)]';
        
        if length(unique(ttail))>2
            tailline=splinefit(ttail,xytail,1,2);
            taillinepts=ppval(tailline,torig(1)-exttail);
            
            xy=[taillinepts xy];
            cptimes=[cptimes(1)-exttail; cptimes];
        end
        thead=cptimes(end-3:end);
        xyhead=[Xgt(thead,id) Ygt(thead,id)]';
        
        if length(unique(thead))>2
            headline=splinefit(thead,xyhead,1,2);
            headlinepts=ppval(headline,torig(end)+exthead);
            
            xy=[xy headlinepts];
            cptimes=[cptimes; cptimes(end)+exthead];
        end
        
%         cptimes=find(Xgt(:,id));
        % cptimes
        
        %     cps
        newmods=newmods+1;
        %     cptimes'
        sfit=splinefit(cptimes,xy,ncps);
        sfit=adjustSplineStruct(sfit, ff, lf, alldpoints, T, 0, [], [], []);
        
        mhs(newmods)=sfit;

    end
end

end