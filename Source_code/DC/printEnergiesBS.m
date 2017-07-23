function energy=printEnergiesBS( i, stateVec, splines,alldpoints, alldvpoints, opt, sceneInfo, labeling, vlabeling, oldEnergy)


global dcStartTime globiter

if opt.verbosity>=3
    [~, ~, energy] = EconBS(stateVec, splines,alldpoints, alldvpoints, opt, sceneInfo, labeling, vlabeling, oldEnergy);

    if i==0
        printMessage(3,'\n  it| time|*to|*ac|*ad|*rm||  Energy|   Data|*Smth|*DetExc| TrjExc|  Lcost| *hreg|  hlin|  hang|  hper|  hocc|   seg||\n');
    end

[totEn, D, S, E, P, L, hreg, hlin, hang, hper, hocc, hseg]= getEnergyValues(energy);

    N=length(splines);


    printMessage(2,'c%3i|%5.1f|  -|%3i|  0|  0||%8.1f|%7.1f|%5.1f|%7.1f|%7.1f|%7.1f|%6.1f|%6.1f|%6.1f|%6.1f|%6.1f|%6.1f|\n', ...
    i, toc(dcStartTime)/60,N, ...
    totEn,D,S,E,P, L,hreg,hlin,hang,hper,hocc,hseg); %%% iter output
    
%     LOG_allens(globiter,:)=[EdetValue EdynValue EexcValue EappValue EperValue EregValue,EoriValue];
end
        
end