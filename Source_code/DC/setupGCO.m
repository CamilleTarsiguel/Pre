function h=setupGCO(nPoints,nLabels,Dcost,Lcost,Scost,Neighborhood)
% set up and return handle for a GCO structure
% to be optimized via GCO alpha expansion
% 



h = GCO_Create(nPoints,nLabels);
GCO_SetDataCost(h,int32(Dcost));
GCO_SetLabelCost(h,int32(Lcost));

if ~isempty(Scost) && any(any(Scost))
    GCO_SetSmoothCost(h,int32(Scost));
    GCO_SetNeighbors(h,Neighborhood);
end

end