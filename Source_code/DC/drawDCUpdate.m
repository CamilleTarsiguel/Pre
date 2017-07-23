function    drawDCUpdate(mhs,used,alldpoints,labeling,outlierLabel,TNeighbors,frames)
% plot the minimization iterations
% 


global opt
if opt.visOptim
    prepFigure; drawPoints(alldpoints,labeling,outlierLabel,TNeighbors);
    drawSplines(mhs,used,labeling,alldpoints,frames)
    drawnow
end


end