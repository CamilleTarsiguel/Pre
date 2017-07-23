function mhs=getSplinesFromEKF(solfile,frames,alldpoints,T)
% Fit splines through an existing EKF solution
% 


mhs=getEmptyModelStruct;
% global tmpcnt processnr

global opt
global tmpcnt processnr
[~,hname]=system('hostname'); hname=hname(1:end-1);

newmods=0;
% for experiment=exp
    % ncps=3;
%     solfile=getSolutionFile(scenario,experiment,'ekf');
%     solfile=fullfile(getHomeFolder,'diss','ekftracking','output',sprintf('s%04d',scenario),sprintf('e%04d.mat',experiment));
    if ~exist(solfile,'file')
        mhs=[];
        
        return;
    end
    
    load(solfile);
%     X=X/1000;Y=Y/1000;
    ncps=max(1,round(length(frameNums)*opt.ncpsPerFrame));
    
    % X=X(1:Fd,:);Y=Y(1:Fd,:);
    % X=removeShorties(X,2);Y=removeShorties(Y,2);
    [F N]=size(X);

    frtokeep=intersect(1:F,frames)';
%     if length(frames)<F && frames(end)<=F
        X=X(frtokeep,:); Y=Y(frtokeep,:);
%     end

    
    
    exthead=2;
    exttail=5;
    [F N]=size(X);
    for id=1:N
        
        
        cptimes=find(X(:,id));
        if length(cptimes)>3
            %% extend to first and last frame
            ff=cptimes(1);  lf=cptimes(end); trl=lf-ff;
            ncps=max(1,round(trl*opt.ncpsPerFrame));
            
            torig=cptimes';
            xy=[X(cptimes,id) Y(cptimes,id)]';
            
            
            ttail=cptimes(1:4);
            xytail=[X(ttail,id) Y(ttail,id)]';
            
            if length(unique(ttail))>2
                tailline=splinefit(ttail,xytail,1,2);
                tailtime=torig(1)-exttail:torig(1)-1;
                taillinepts=ppval(tailline,tailtime);
                
                xy=[taillinepts xy];
                cptimes=[(cptimes(1)-exttail:cptimes(1)-1)'; cptimes];
            end
            thead=cptimes(end-3:end);
            xyhead=[X(thead,id) Y(thead,id)]';
            
            if length(unique(thead))>2
                headline=splinefit(thead,xyhead,1,2);
                headlinepts=ppval(headline,torig(end)+1:torig(end)+exthead);
                
                xy=[xy headlinepts];
                cptimes=[cptimes; (cptimes(end)+1:cptimes(end)+exthead)'];
            end
            
%              if tmpcnt<2000, save(sprintf('tmp/tmp%d/tmpvar_%05d.mat',processnr,tmpcnt),'*'); tmpcnt=tmpcnt+1; end
            
            newmods=newmods+1;
            sfit=splinefit(cptimes,xy,ncps);
            sfit=adjustSplineStruct(sfit,ff,lf,alldpoints, F, 0, [], [], []);
%             sfit.start=ff; sfit.end=lf;
%             [sfit.labelCost, sfit.lcComponents]=getSplineGoodness(sfit,1,alldpoints,F);
%             sfit.lastused=0;
            mhs(newmods)=sfit;
        end
% 		save(sprintf('tmp/tmp%d/tmpvar_%05d.mat',processnr,tmpcnt),'*'); tmpcnt=tmpcnt+1;

    end
% end

%% check for doubles

if ~isempty(mhs)
    dists=Inf*ones(newmods);
    for m1=1:newmods
        for m2=m1+1:newmods
            if mhs(m1).pieces==mhs(m2).pieces
                dists(m1,m2)=sum(sum(abs(mhs(m1).coefs-mhs(m2).coefs)));
            end
        end
    end
    
    redundant=sum(~dists,2);
    mhs=mhs(~redundant);
end
% save(sprintf('tmp/tmp%d/tmpvar_%05d.mat',processnr,tmpcnt),'*'); tmpcnt=tmpcnt+1;

end