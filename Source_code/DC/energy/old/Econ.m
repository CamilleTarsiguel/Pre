function [fx, dfx, energy]= Econ( ...
    coefs,splines,alldpoints, alldvpoints, opt, sceneInfo, labeling, vlabeling, oldEnergy)

% This function computes the energy of the current state
% and the derivatives w.r.t. the spline coefficients
%
% Returns the energy value fx and its derivative dfx
% energy is a struct containing individual components

global stateVec parEper spInfo br indexexc ptindex gspl index
gspl=splines;

stateVec=coefs;

% global EdatValue ElinValue EangValue EsegValue
global ptsv pts ptInfo

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
    
% mex stuff
if opt.mex
    spInfo=[[splines(:).pieces];[splines(:).start];[splines(:).end]]';

    br=[splines(:).breaks];
    index=[splines(:).index];

    ptindex=[splines(:).ptindex];
    
    brInfo=ones(length(splines),2);
    for id=1:length(splines)
        brInfo(id,1)=splines(id).breaks(1);
        brInfo(id,2)=splines(id).breaks(end);
    end
        
    pcindex=[];
    brindex=[]; brindxcnt=1;
    Nsp=length(splines);
    for id=1:Nsp
        pcindex=[pcindex id*ones(1,splines(id).pieces)];
        brindex=[brindex brindxcnt:brindxcnt+splines(id).pieces-1];
        brindxcnt=brindxcnt+splines(id).pieces+1;
    end
    indexexc=[pcindex; brindex];
    indexexc=indexexc-1;

end

%%%% FOR TESTING BSPLINES
%  spl = vecToSplines(stateVec,splines);
%  
%  
%  for id=1:length(splines)
%      bspl(id)=mypp2sp(spl(id));
%  %     spl(id)
%  %     bspl(id)
%  end
%  bsv=BSplinesToVec(bspl);

%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edat  -  data term  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EdatValue=0; dEdat=zeros(length(stateVec),1);
if wtEdat
    if opt.mex
        [EdatValue, dEdat]=Edat_mex(stateVec,parEdat,spInfo,br,ptindex,ptsv);
    else
        if nargout>1
            [EdatValue, dEdat]=Edat(stateVec,splines,parEdat, alldvpoints, vlabeling);
        else
            EdatValue=Edat(stateVec,splines,parEdat, alldvpoints, vlabeling);
        end
    end
end
%  [EdatValueBS, dEdatBS]=EdatBS(bsv,bspl,parEdat, alldvpoints, vlabeling);
%  [EdatValue, EdatValueBS]

% EdatValue
% pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elin  -  linear velocity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global glsplines
% glsplines=splines;
ElinValue=0; dElin=zeros(length(stateVec),1);
if wtElin
    if opt.mex
        [ElinValue, dElin]=Elin_mex(stateVec,parElin,spInfo,br,index);
    else
        if nargout>1
            [ElinValue, dElin]=Elin(stateVec,splines,parElin);
        else
            ElinValue=Elin(stateVec,splines,parElin);
        end
    end
end

%  [ElinValueBS, dElinBS]=ElinBS(bsv,bspl,parElin);
%  [ElinValue, ElinValueBS]

% pause
% global stateVec parElin spInfo br index splines
% pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eang  -  angular velocity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EangValue=0; dEang=zeros(length(stateVec),1);
if wtEang
     if opt.mex
         [EangValue, dEang]=Eang_mex(stateVec,parEang,spInfo,br,index);
     else
        if nargout>1
            [EangValue, dEang]=Eang(stateVec,splines);
        else
            EangValue=Eang(stateVec,splines);
        end
     end
end
%  [EangValueBS, dEangBS]=EangBS(bsv,bspl);
%  [EangValue, EangValueBS]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eper  -  persistence   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EperValue=0; dEper=zeros(length(stateVec),1);
if wtEper
    if opt.mex
        [EperValue, dEper]=Eper_mex(stateVec,parEper,spInfo,br,index,sceneInfo.imOnGP);
    else
        if nargout>1
            [EperValue, dEper]=Eper(stateVec,splines, parEper);
        else
            EperValue=Eper(stateVec,splines, parEper);
        end
    end
end
%  [EperValueBS, dEperBS]=EperBS(bsv,bspl, parEper);
%  [EperValue, EperValueBS]
% pause
% [EperValuemex, dEpermex]=Eper_mex(stateVec,parEper,spInfo,br,index,sceneInfo.imOnGP);
% [EperValuemat, dEpermat]=Eper(stateVec,splines, parEper);
% assert(abs(EperValuemex-EperValuemat)<1e-5,'%f %f',EperValuemex,EperValuemat);
% assert(sum(abs(dEpermex-dEpermat))<1e-2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eexc -  traj.-level exclusion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EexcValue=0; dEexc=zeros(length(stateVec),1);
if wtEexc
    if opt.mex
        [EexcValue, dEexc]=Eexc_mex(stateVec,parEexc,spInfo,brInfo,br,indexexc);
    else
        if nargout>1
            [EexcValue, dEexc]=Eexc(stateVec,splines);
        else
            EexcValue=Eexc(stateVec,splines);
        end
    end
end
%  [EexcValueBS, dEexcBS]=EexcBS(bsv,bspl);
%  [EexcValue, EexcValueBS]

% save('tmpexcvar.mat','stateVec','parEexc','spInfo','brInfo','br','indexexc');
% pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Efid -  data fidelity  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EfidValue=0; dEfid=zeros(length(stateVec),1);
if wtEfid
    if opt.mex
        [EfidValue, dEfid]=Efid_mex(stateVec,parEfid,spInfo,br,index,pts, ptInfo);
    else
        if nargout>1
            [EfidValue, dEfid]=Efid(stateVec,splines, parEfid, alldpoints);
        else
            EfidValue=Efid(stateVec,splines, parEfid,alldpoints);
        end
    end
end
% pause
% [EfidValuemex, dEfidmex]=Efid_mex(stateVec,parEfid,spInfo,br,index,pts,ptInfo);
% [EfidValuemat, dEfidmat]=Efid(stateVec,splines, parEfid);
% assert(abs(EfidValuemex-EfidValuemat)<1e-5,'%f %f',EfidValuemex,EfidValuemat);
% assert(sum(abs(dEfidmex-dEfidmat))<1e-2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eseg  -  segments consistency %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EsegValue=0; dEseg=zeros(length(stateVec),1);
if wtEseg
    if opt.mex
        [EsegValue, dEseg]=Eseg_mex(stateVec,parEseg,spInfo,br,index);
    else
        if nargout>1
            [EsegValue, dEseg]=Eseg(stateVec,splines,opt);
        else
            EsegValue=Eseg(stateVec,splines,opt);
        end
    end
end
%  pause(.1);

% apply weighting
EdatValue = wtEdat * EdatValue; dEdat = wtEdat * dEdat;
ElinValue = wtElin * ElinValue; dElin = wtElin * dElin;
EangValue = wtEang * EangValue; dEang = wtEang * dEang;
EperValue = wtEper * EperValue; dEper = wtEper * dEper;
EexcValue = wtEexc * EexcValue; dEexc = wtEexc * dEexc;
EfidValue = wtEfid * EfidValue; dEfid = wtEfid * dEfid;
EsegValue = wtEseg * EsegValue; dEseg = wtEseg * dEseg;




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

% detExc is constant
% ?????????



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