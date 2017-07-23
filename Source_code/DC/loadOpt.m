function opt=loadOpt(opt,fname)
% load options from file
allopts=load(fname);

opt.outlierCost=allopts(1);
opt.labelCost=allopts(2);
opt.unaryFactor=allopts(3);
opt.pairwiseFactor=allopts(4);
opt.goodnessFactor=allopts(5);
opt.fidelityFactor = allopts(6);
opt.persistenceFactor = allopts(7);
opt.lengthFactor = allopts(8);
opt.exclusionFactor = 0;
opt.proxcostFactor = 0;
opt.curvatureFactor = allopts(11);
opt.slopeFactor = allopts(12);

end
