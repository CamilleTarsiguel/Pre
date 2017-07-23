function [metrics2d, metrics3d]= ...
    printDCUpdateOGM(stateInfo,splines,used,nNew,nRemoved,itcnt,allEnergies,itType,titlestr)

% print energy and performance information
% 


global opt gtInfo dcStartTime sceneInfo

gtheader='';
metheader='';
metrics2d=zeros(1,14);
metrics3d=zeros(1,14);
% if nargin<10
%     itType='';
% end

if sceneInfo.gtAvailable && opt.verbosity>=3 
    gtheader=''; metheader='';
    
    if opt.track3d
        gtheader='  ----------- M E T R I C S (3D)---------- |||';
        metheader=' MOTA  MOTP| GT  MT  ML|  FP   FN IDs  FM  |||';
        if opt.met2d
            gtheader=[gtheader '  ----------- M E T R I C S (2D)--------- '];
            metheader=[metheader ' MOTA  MOTP| GT  MT  ML|  FP   FN IDs  FM|'];
        end
    else
        if opt.met2d
            gtheader='  ----------- M E T R I C S (2D)---------- |||';
            metheader=' MOTA  MOTP| GT  MT  ML|  FP   FN IDs  FM |';
        end        
    end
    
end

if ~mod(itcnt,10)
    printMessage(2,'\n -- S: %04d, (%3d : %3d) --|| ------------- ENERGY  VALUES -----------------| ---------------- Label Cost ------------ |||%s',sceneInfo.scenario,sceneInfo.frameNums(1), sceneInfo.frameNums(end),gtheader);
    printMessage(2,'\n   it| time|tot|act|add|rem||  Energy|    Data| Smth| DetExc| TrjExc|  Lcost|   hreg|  hlin|  hang|  hper|  hocc|  hseg|||%s\n',metheader);
end

[totEn, D, S, E, P, L, hreg, hlin, hang, hper, hocc, hseg] = ...
    getEnergyValues(allEnergies);

%                 tpcnt| time|tot|act|add|rem|| Ener|  Data| Smth| DetE| TrjE| Lcst|  hreg| hlin| hang| hper| hocc| hseg|
printMessage(2,' %1s%3i|%5.1f|%3i|%3i|%3i|%3i||%8.1f|%8.1f|%5.1f|%7.1f|%7.1f|%7.1f|%7.1f|%6.1f|%6.1f|%6.1f|%6.1f|%6.1f|||', ...
    itType,itcnt, toc(dcStartTime)/60,length(splines),length(used),nNew, nRemoved, ...
    totEn,D,S,E,P,L,hreg,hlin,hang,hper,hocc,hseg); %%% iter output
% pause

if opt.verbosity>=3
            
    if sceneInfo.gtAvailable
    
	% get state (x,y) from active splines
        stateInfo=getStateFromSplines(splines(used), stateInfo);
        
        % cut state to tracking area if needed
        if opt.cutToTA
            stateInfo=cutStateToTrackingArea(stateInfo);
        end

        % compute 3d metrics
        if opt.track3d
            [metrics3d, metrNames3d]=CLEAR_MOT(gtInfo,stateInfo,struct('eval3d','1'));
            printMetrics(metrics3d,metrNames3d,0,[12 13 4 5 7 8 9 10 11]);
            printMessage(3,'|||');
        end    

        % compute 2d metrics
        if opt.met2d && (~mod(itcnt,5) || itcnt==1) % only display 2d metrics every 5 iterations
            if opt.track3d
                [stateInfo.Xi, stateInfo.Yi]=projectToImage(stateInfo.X, stateInfo.Y,sceneInfo);
            else
                stateInfo.Xi=stateInfo.X; stateInfo.Yi=stateInfo.Y;
            end
            stateInfo=getBBoxesFromState(stateInfo);
            
            % YSHIFT
            if sceneInfo.yshift
                stateInfo.Yi=stateInfo.Yi+stateInfo.H/2;
            end
            
%             evopt.eval3d=0;
%             [metrics2d, metricsInfo2d]=CLEAR_MOT(gtInfo,stateInfo,evopt);
%             printMetrics(metrics2d,metricsInfo2d,0,[12 13 4 5 7 8 9 10 11]);    
        end

    end
end
printMessage(2,'\n');


if exist('titlestr','var') && opt.visOptim
    mota=metrics2d(12);
    if opt.track3d,        mota=metrics3d(12); end
    title(sprintf('E: %d, MOTA: %.2f%%',int32(D)+int32(S)+int32(L)+int32(E),mota));
    
end

end
