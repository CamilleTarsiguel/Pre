function energy = energyLogs(energy, itType, itCnt, tocIt, tocGlobal)
% insert more info into energy struct

energy.itType=itType;
energy.itCnt=itCnt;
energy.localTime=tocIt;
energy.globalTime=tocGlobal;
end