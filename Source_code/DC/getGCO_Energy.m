function [inE inD inS inL] = getGCO_Energy(Dcost, Scost, Lcost, TNeighbors, labeling)
% compute energy values for a given labeling (Graph Cuts Optimization GCO)
% 


[nLabels nPoints]=size(Dcost);
h = GCO_Create(nPoints,nLabels);
GCO_SetDataCost(h,int32(Dcost));
GCO_SetLabelCost(h,int32(Lcost));
if ~isempty(Scost)
    GCO_SetSmoothCost(h,int32(Scost));
    GCO_SetNeighbors(h,TNeighbors);
end

GCO_SetLabeling(h,int32(labeling));

[inE inD inS inL] = GCO_ComputeEnergy(h);
end
