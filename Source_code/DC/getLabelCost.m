function Lcost=getLabelCost(splines,opt)
% Compute the label cost for all splines
% The label cost is not computed but just copied here
% It is mainly computed in getSplineGoodness.m
% 


% global opt

Lcost=0;
if isempty(splines)
    return;
end

nCurModels=length(splines);

lc1=opt.labelCost*ones(1,nCurModels); % standard uniform
lc2=[splines(:).labelCost]; % hdyn, hfid, hper,...
lc3=0;%proxcostFactor*proxcost;
Lcost =[lc1+lc2+lc3 0]; % extra zero for the outlier, it has no label cost
    
end