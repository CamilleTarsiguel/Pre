%% used in supplemental video

ensnapsdir='ensnaps/visoptim';
allfiles=dir(fullfile(ensnapsdir,'*.mat'));
nfiles=length(allfiles);

envals=zeros(0,2);
enCnt=0;
tictocs=zeros(size(envals));
mets2d=zeros(0,2,14);
mets3d=zeros(0,2,14);

dotSize=10;
lw=2;
global sceneInfo opt

addpath(genpath('../../contracking/utils/'));
addpath('mex/bin/'); 
selectedfiles=1:nfiles;
selectedfiles=2670:3391;
selectedfiles=9144:9688;
selectedfiles=10099:10539;
selectedfiles=10360:10539;
for n=selectedfiles
    n
    fname=fullfile(ensnapsdir,allfiles(n).name)
    load(fname);
    labeling(labeling==-1)=length(splines)+1;
%     length(splines)
%     unique(labeling)
%     pause
%     energy=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);
    
    F=opt.seqLength;
    stateInfo.F=F;
    stateInfo=getStateFromSplines(splines,stateInfo,1);
    X=stateInfo.Xi; Y=stateInfo.Yi;
%     X(1:2,:)
%     pause
    [~, N]=size(X);
    % how many splines?
    nCurModels=length(splines);
    nLabels=nCurModels+1;
    outlierLabel=nLabels;
    used=setdiff(unique(labeling),outlierLabel);
    
    newlabeling=labeling;
    newlabeling(labeling==outlierLabel)=-1;
    
    detections=setDetectionsIDs(detections,newlabeling);

    % display static image
    t=round(opt.seqLength/2);
    clf
    im=imread([sceneInfo.imgFolder sprintf(sceneInfo.imgFileFormat,sceneInfo.frameNums(t))]);
    if (size(im,3)==1), im=repmat(im,[1 1 3]); end % greyscale
%     im=ones(sceneInfo.imgHeight, sceneInfo.imgWidth); %
    imshow(im,'Border','tight')
    hold on
%     pause

    %%% show dots
    for t=1:F
        nboxes=length(detections(t).xp);        
        
        
        % foot position (ordering: outliers, actives)
        outliers=find(detections(t).id==-1);
        detcol=ones(1,3);
        plot(detections(t).xi(outliers),detections(t).yi(outliers), ...
                     'o','color',detcol,'MarkerSize',10);
%         for k=1:nboxes
%             if detections(t).id(k) == -1
% %                 pause
%                 detcol=ones(1,3);
%                 plot(detections(t).xi(k),detections(t).yi(k), ...
%                     'o','color',detcol,'MarkerSize',dotSize*max(0.1,detections(t).sc(k)));
%             end
%         end
        for k=1:nboxes
            if detections(t).id(k) ~= -1
                detcol=getColorFromID(detections(t).id(k));
                plot(detections(t).xi(k),detections(t).yi(k), ...
                    '.','color',detcol,'MarkerSize',dotSize*max(0.1,detections(t).sc(k)));
            end
        end
    end
    
    %%% show trajectories (ordering: outliers, actives)
    for id=1:N
        extar=find(X(:,id));
        % inactives are white
        if isempty(intersect(used,id))
            col=ones(1,3);
        
            plot(X(extar,id),Y(extar,id),'color',col, ...
                'linewidth',lw/2);
        end
        
    end
%     pause
    for id=1:N
        extar=find(X(:,id));
        col=getColorFromID(id);
        if ~isempty(intersect(used,id))
            plot(X(extar,id),Y(extar,id),'color',col, ...
                'linewidth',lw);

        end
    end
%     pause
    
    %%% text energy
    % get energy
    energy=evaluateEnergy(alldpoints, Nhood, labeling, splines, opt, []);
%     text(20,500,sprintf('Energy:%.1f',energy.value),'FontSize',20);
%     text(20,540,sprintf('Iteration type: %d', itType),'FontSize',20);
    im2save=getframe(gcf);
    im2save=im2save.cdata;
    imwrite(im2save,[fname '.png']);
    pause(.01);
end