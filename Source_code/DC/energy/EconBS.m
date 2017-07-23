% function [fx, dfx, energy]= EconBS( ...
%     coefs,splines,alldpoints, alldvpoints, opt, sceneInfo, labeling, vlabeling, oldEnergy)

function [fx, dfx, energy, fxDat, dfxDat, fxLin, dfxLin, fxAng, dfxAng, ...
    fxPer, dfxPer, fxExc, dfxExc, fxSeg, dfxSeg]= EconBS( ...
    coefs,splines,alldpoints, alldvpoints, opt, sceneInfo, labeling, vlabeling, oldEnergy)

% This function computes the energy of the current state
% and the derivatives w.r.t. the spline coefficients
%
% Returns the energy value fx and its derivative dfx
% energy is a struct containing individual components


fx=0;
dfx=zeros(length(coefs),1);
energy=oldEnergy;
stateVec=coefs;

wtEdat=opt.unaryFactor;
wtElin=opt.slopeFactor;
wtEang=opt.curvatureFactor;
wtEper=opt.persistenceFactor;
wtEexc=opt.proxcostFactor;
wtEfid=opt.fidelityFactor;
wtEseg=opt.segFactor;

% wtEdat=0; wtElin=1; wtEang=0; wtEseg=0; 

% parameters
parEdat=opt.conOpt.enParEdat; % Edat par
parElin=opt.conOpt.enParElin; % Elin par: normfac, m/s, fR
parEang=opt.conOpt.enParEang; % Eang par
parEper=opt.conOpt.enParEper; % Eper par
parEexc=opt.conOpt.enParEexc; % Eexc par
parEfid=opt.conOpt.enParEfid; % Efid par
parEseg=opt.conOpt.enParEseg; % Eseg par

% save bounding boxes
bm=sceneInfo.targetSize*2;
tmpspl=vecToBSplines(coefs,splines);
for id=1:length(splines)
    tsp=tmpspl(id);
    xyvals=spval(tsp,tsp.start:tsp.end);
    splines(id).bbox=[min(xyvals(1,:))-bm max(xyvals(1,:))+bm min(xyvals(2,:))-bm max(xyvals(2,:))+bm];
end

% MEX STUFF
% offsets: where in the state vector does each spline start?? 0-based
sna=[splines(:).number]*2;
sna=cumsum(sna)-sna;

spInfo=[sna;[splines(:).start];[splines(:).end];[splines(:).number]-3]';
knots=[splines(:).knots];
% br=[splines(:).breaks];
index=[splines(:).index];
ptindex=[splines(:).ptindex];


global ptsv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edat  -  data term  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EdatValue=0; dEdat=zeros(length(stateVec),1);
if wtEdat
    if opt.mex      
%         save('EdatBS.mat','stateVec','parEdat','spInfo','knots','ptindex','ptsv','splines','alldvpoints','vlabeling');
%         fprintf('a');
        [EdatValue, dEdat]=EdatBS_mex(stateVec,parEdat,spInfo,knots,ptindex,ptsv);
%         fprintf('b');
    else
        if nargout>1
            [EdatValue, dEdat]=EdatBS(stateVec,splines,parEdat, alldvpoints, vlabeling);            
        else
            EdatValue=EdatBS(stateVec,splines,parEdat, alldvpoints, vlabeling);
        end
    end
end
fxDat=EdatValue; dfxDat=dEdat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elin  -  linear velocity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global glsplines
% glsplines=splines;
ElinValue=0; dElin=zeros(length(stateVec),1);
if wtElin
    if opt.mex
%         save('ElinBS.mat','stateVec','parElin','spInfo','knots','index');
        [ElinValue, dElin]=ElinBS_mex(stateVec,parElin,spInfo,knots,index);
    else
        if nargout>1
            [ElinValue, dElin]=ElinBS(stateVec,splines,parElin);            
        else
            ElinValue=ElinBS(stateVec,splines,parElin);
            
        end
    end
end
fxLin=ElinValue; dfxLin=dElin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eang  -  angular velocity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EangValue=0; dEang=zeros(length(stateVec),1);
if wtEang
     if opt.mex
%          save('EangBS.mat','stateVec','parEang','spInfo','knots','index');
         [EangValue, dEang]=EangBS_mex(stateVec,parEang,spInfo,knots,index);
%          [EangValue, dEang]=EangBS(stateVec,splines,parEang);            
     else
        if nargout>1
            [EangValue, dEang]=EangBS(stateVec,splines,parEang);            
        else
            EangValue=EangBS(stateVec,splines,parEang);
        end
     end
end
fxAng=EangValue; dfxAng=dEang;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eper  -  persistence   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EperValue=0; dEper=zeros(length(stateVec),1);
if wtEper
    if opt.mex        
%         save('EperBS.mat','stateVec','parEper','spInfo','knots','index','sceneInfo');
        [EperValue, dEper]=EperBS_mex(stateVec,parEper,spInfo,knots,index,sceneInfo.imOnGP);
    else
        if nargout>1
            [EperValue, dEper]=EperBS(stateVec,splines, parEper);            
        else
            EperValue=EperBS(stateVec,splines, parEper);
        end
    end
end
fxPer=EperValue; dfxPer=dEper;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eexc -  traj.-level exclusion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EexcValue=0; dEexc=zeros(length(stateVec),1);
if wtEexc
    if 0%opt.mex
        %TODO
        [EexcValue, dEexc]=EexcBS_mex(stateVec,parEexc,spInfo,brInfo,knots,indexexc);
    else
        if nargout>1
            [EexcValue, dEexc]=EexcBS(stateVec,splines, parEexc);
        else
            EexcValue=EexcBS(stateVec,splines,parEexc);
        end
    end
end
fxExc=EexcValue; dfxExc=dEexc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Efid -  data fidelity  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EfidValue=0; dEfid=zeros(length(stateVec),1);
% this term, originally used in CVPR 2012 and CVPR 2013 papers,
% has been dropped for the PAMI publication

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eseg  -  segments consistency %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EsegValue=0; dEseg=zeros(length(stateVec),1);
if wtEseg
    if opt.mex
%          save('EbndBS.mat','stateVec','parEseg','spInfo','knots','index');
        [EsegValue, dEseg]=EbndBS_mex(stateVec,parEseg,spInfo,knots,index);
    else
        if nargout>1
            [EsegValue, dEseg]=EbndBS(stateVec,splines,opt);
        else
            EsegValue=EbndBS(stateVec,splines,opt);
        end
    end
end
fxSeg=EsegValue; dfxSeg=dEseg;

% apply weighting
EdatValue = wtEdat * EdatValue; dEdat = wtEdat * dEdat;
ElinValue = wtElin * ElinValue; dElin = wtElin * dElin;
EangValue = wtEang * EangValue; dEang = wtEang * dEang;
EperValue = wtEper * EperValue; dEper = wtEper * dEper;
EexcValue = wtEexc * EexcValue; dEexc = wtEexc * dEexc;
EfidValue = wtEfid * EfidValue; dEfid = wtEfid * dEfid;
EsegValue = wtEseg * EsegValue; dEseg = wtEseg * dEseg;

% the following energy components remain constant during
% continuous optimization

% outlier cost is constant
outlierLabel=-1;
outliers=find(labeling==outlierLabel);
dataOutliers = opt.outlierCost*sum(alldpoints.sp(outliers));
EdatValue = EdatValue + dataOutliers;

% label cost is constant
ElabValue = opt.labelCost * length(splines);

% global oldEnergy
% Smoothness is constant
EsmtValue = oldEnergy.smoothness;
EdetexcValue = oldEnergy.detExclusion;



% final value is a linear combination of all terms
fx = ... 
    EdatValue + ...
    ElinValue + ...
    EangValue + ...
    EperValue + ...
    EexcValue + ...
    EfidValue + ...
    EsegValue + ... % rest is constant during continuous optim
    ElabValue + ...
    EsmtValue;


% and the gradient
if nargout>1
    dfx = ...
        dEdat + ...
        dElin + ...
        dEang + ...
        dEper + ...
        dEexc + ...
        dEfid + ...
        dEseg;    

end

% set struct
energy=oldEnergy;
if nargout>2
    energy = setEnergyValues(EdatValue, EsmtValue, EdetexcValue, ...
        EexcValue, ElabValue, ElinValue, EangValue, EperValue, EfidValue, EsegValue);    
end


end