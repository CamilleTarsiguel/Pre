function [totEn, D, S, E, P, L, hreg, hlin, hang, hper, hocc, hseg]= ...
    getEnergyValues(energy)
% retrieve individual energy components from struct

totEn = energy.value;
D = energy.data;
S = energy.smoothness;
E = energy.detExclusion;
P = energy.trajExclusion;
L = energy.labelCost.value;
hreg = energy.labelCost.regularizer;
hlin = energy.labelCost.linVelocity;
hang = energy.labelCost.angVelocity;
hper = energy.labelCost.persistence;
hocc = energy.labelCost.occlusionGaps;
hseg = energy.labelCost.segConsistency;

enTol=1e-5;

assert(abs(totEn - (D+S+E+P+L)) < enTol, 'energy value is wrong');
% energy
% energy.labelCost
% (hreg+hlin+hang+hper+hocc)
assert(abs(L-(hreg+hlin+hang+hper+hocc+hseg)) < enTol,'label cost value is wrong');

end