function [E, D, S, L, labeling] = ...
    doAlphaExpansion(Dcost, Scost, Lcost, Neighborhood)

% minimize E(f,T) wrt. f by alpha expansion
% Note: This only works for submodular energies
% 
% The gco code is available at
% http://vision.csd.uwo.ca/code/

[nLabels, nPoints]=size(Dcost);

% trivial case for one label
if nLabels==1
    D=sum(Dcost);    S=0; L=Lcost;    E=D+S+L;
    labeling=ones(1,nPoints);
    return;
end


% set up GCO structure

h=setupGCO(nPoints,nLabels,Dcost,Lcost,Scost,Neighborhood);

GCO_SetLabelOrder(h,1:nLabels);
GCO_Expansion(h);
labeling=GCO_GetLabeling(h)';
[E, D, S, L] = GCO_ComputeEnergy(h);

% clean up
GCO_Delete(h);
end