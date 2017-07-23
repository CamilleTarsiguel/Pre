function printHeader(sceneInfo,scenario,randrun)


printMessage(2,' =====================================================================\n');
printMessage(2,'|                Discrete-Continuous Optimization                     |\n');
printMessage(2,'|                With Explicit Exclusion Handling                     |\n');
printMessage(2,'|                                                                     |\n');
if isnumeric(scenario)
    printMessage(2,'|       Scenario: %10d              Random Run: %12d    |\n', ...
        scenario, randrun);
elseif ischar(scenario)
    printMessage(2,'|  Scenario: %20s      Random Run: %15d    |\n', ...
        scenario, randrun);
end

if all(isfield(sceneInfo,{'dataset','sequence'}))
printMessage(2,'|  Dataset: %21s           Sequence: %12s    |\n',sceneInfo.dataset,sceneInfo.sequence);
end

% :\nSCENARIO %i, LEXP %i, RAND %i: %i frames (%i:%i:%i)\n', ...
%     scenario,lexperiment,randrun,length(frameNums),frameNums(1),frameNums(2)-frameNums(1),frameNums(end));
printMessage(2,' =====================================================================\n\n');

end