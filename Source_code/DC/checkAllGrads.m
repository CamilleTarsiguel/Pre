%%checkAllGrads
splines=splinesActive;

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

%%
% d=checkgrad('Econ',xvec,1e-3,splinesActive,ptsAndVPoints,dcOpt,sceneInfo,newLabeling,energy)
% d=checkgrad('Edat',xvec,1e-3,splinesActive,dcOpt.conOpt.enParEdat,alldpoints,newLabeling)
d=checkgrad('Edat_mex',xvec,1e-3,dcOpt.conOpt.enParEdat,spInfo,br,ptindex,ptsv)
% d=checkgrad('Elin',xvec,1e-3,splinesActive,dcOpt.conOpt.enParElin)
d=checkgrad('Elin_mex',xvec,1e-3,dcOpt.conOpt.enParElin,spInfo,br,index)
% d=checkgrad('Eang',xvec,1e-3,splinesActive)
d=checkgrad('Eang_mex',xvec,1e-3,dcOpt.conOpt.enParElin,spInfo,br,index)
% d=checkgrad('Eper',xvec,1e-3,splinesActive,dcOpt.conOpt.enParEper)
d=checkgrad('Eper_mex',xvec,1e-3,dcOpt.conOpt.enParEper,spInfo,br,index,sceneInfo.imOnGP)
% d=checkgrad('Eexc',xvec,1e-3,splinesActive)
d=checkgrad('Eexc_mex',xvec,1e-3,dcOpt.conOpt.enParEexc,spInfo,brInfo,br,indexexc)
% d=checkgrad('Efid',xvec,1e-3,splinesActive,dcOpt.conOpt.enParEfid,alldpoints)
d=checkgrad('Efid_mex',xvec,1e-3,dcOpt.conOpt.enParEfid,spInfo,br,index,pts,ptInfo)