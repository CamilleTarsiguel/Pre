function energy = setEnergyValues(Edat, Esmt, Edetexc, ...
    Etrjexc, hreg, hlin, hang, hper, hocc, hseg)

% form a struct containing all energy components

labcost= hreg + hlin + hang + hper + hocc + hseg;
totEn = Edat + Esmt + Edetexc + Etrjexc + labcost;

energy.value = totEn;
energy.data = Edat;
energy.smoothness = Esmt;
energy.detExclusion = Edetexc;
energy.trajExclusion = Etrjexc;
energy.labelCost.value = labcost;
energy.labelCost.regularizer = hreg;
energy.labelCost.linVelocity = hlin;
energy.labelCost.angVelocity = hang;
energy.labelCost.persistence = hper;
energy.labelCost.occlusionGaps = hocc;
energy.labelCost.segConsistency = hseg;