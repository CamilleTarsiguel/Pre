function opt=getAuxOpt(conffile,opt,sceneInfo,T)
% fill auxiliary options
%  These mainly control parameters for the continuous 
%  approximations of the labelcost, frame rate, etc...

if exist(conffile,'file')
    auxOpt=load(conffile);
else
    auxOpt=[0.1 1 0.05 1 1 2 0 0];
end

nf=1; if opt.track3d, nf=1000; end	% normalization factor, mm for 3d tracking

if sceneInfo.scenario == 51 || sceneInfo.scenario == 53, auxOpt(7)=.1; auxOpt(8)=250; end

opt.conOpt.enParEdat=[opt.dataFunction, nf, auxOpt(1)]; % Edat par
opt.conOpt.enParElin=[nf, auxOpt(2), sceneInfo.frameRate, auxOpt(7), auxOpt(8)]; % Elin par: normfac, m/s, fR, epsilon, ypos of horizon
opt.conOpt.enParEang=[nf, sceneInfo.frameRate,  auxOpt(7), auxOpt(8)]; % Eang par, ...,  epsilon, ypos of horizon
opt.conOpt.enParEper=[nf, opt.borderMargin/nf, T]; % Eper par
opt.conOpt.enParEexc=[auxOpt(3), sceneInfo.targetSize nf]; % Eexc par
opt.conOpt.enParEfid=[auxOpt(4) sceneInfo.targetSize/2 auxOpt(5)]; % Efid par
opt.conOpt.enParEseg=auxOpt(6); % Delta
opt.seqLength=T;


end