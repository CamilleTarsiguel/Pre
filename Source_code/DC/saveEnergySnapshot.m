function saveEnergySnapshot(alldpoints, detections, Nhood, labeling, ...
    splines, opt, sceneInfo, itType)
% for debugging and evaluation, save everything relevant for each iteration

if ~opt.saveEnSS
	return;
end


    
timestamp=datestr(now,'yyyy-mm-dd_HH-MM-SS-FFF');
fname=['ensnaps/visoptim/' timestamp '-' num2str(itType) '.mat'];

save(fname,'*');
% pause(.5);

end