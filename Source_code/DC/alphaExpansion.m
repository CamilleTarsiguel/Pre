function [Eogm, logm]=alphaExpansion(Dcost,Scost,Lcost, labeling, ...
    Nhood, alldpoints, mhs,   dcOpt)
% alpha expansion with openGM as inference engine
%
TNeighbors=Nhood.TNeighbors;
SNeighbors=Nhood.SNeighbors;

%% some special cases...
% do Olga's code if energy is submodular
if ~dcOpt.exclusionFactor && ~dcOpt.proxcostFactor
    [E, ~, ~, ~, olga_labeling]=doAlphaExpansion(Dcost, Scost, Lcost, TNeighbors);
    logm=double(olga_labeling);    Eogm=double(E);
    return;
end

% if we only have unaries, inference is trivial
if dcOpt.unaryFactor && ~dcOpt.pairwiseFactor && ~dcOpt.exclusionFactor && ~sum(Lcost)
    [m, mm]=min(Dcost);
    logm=mm;
    Eogm=sum(m);
    return;
end

% if we have only one (outlier) label
if size(Dcost,1)<2
    logm=labeling;
    Eogm=length(labeling)*dcOpt.outlierCost;
    return;
end

%% The general case starts here
global DcostAlphaAux

% first determine algorithm
switch(dcOpt.disOpt.alg)
    case 'TRWS'
        disAlg=1;
    case 'QPBO'
        disAlg=2;
    case 'MQPBO'
        disAlg=3;
    case 'ICM'
        disAlg=4;
    case 'TRWSi'
        disAlg=5;
    case 'FastPD'
        disAlg=6;
end
%%


% for debugging...
global Auxmatcl Lcostcl LcostInd binLabeling Dpw Dpwi Expw Expwi  notAlphasNoOutlier pwLcost nAux

% for energy snapshots
global detections sceneInfo
sscnt=0;

T=dcOpt.seqLength;
tmpState.F=T;
tmpState=getStateFromSplines(mhs, tmpState, ~dcOpt.track3d);
[prox, proxt, proxcost]=getSplineProximity(mhs,T,[],dcOpt);
tmpState2=[];

% precompute energy values that remain constant
precomp.Dcost=Dcost;
precomp.Lcost=Lcost;
precomp.proxcost=proxcost;
% get state for all splines over entire sequence
tt=1:T;
for id=1:length(mhs)
    allxy=ppval(mhs(id),tt)';
    tmpState2.X(tt,id)=allxy(:,1);    tmpState2.Y(tt,id)=allxy(:,2);
end

% do not use Lazy Flipper by default
if ~isfield(dcOpt,'useLF'), dcOpt.useLF=0; end

% prox
% proxcost
% tmpproxcost=proxcost+proxcost'
% andlabelcost=~~proxcost'

% what is the LabelSpace?
[nLabels, nPoints]=size(Dcost);
outlierLabel=nLabels;   % the highest label is used for the outlier model


% ordering for alpha expansion
labelOrder=1:nLabels;
% labelOrder=randperm(nLabels); % random seems to perform worse.

% initial labeling, just set all to outliers
% labeling=nLabels*ones(1,nPoints);

% graph cuts of a simplified submodular energy
% if dcOpt.disOpt.initSimple
[~, ~, ~, ~, olga_labeling]=doAlphaExpansion(Dcost, Scost, Lcost, TNeighbors);
% end

% initial labeling energy
Pw=TNeighbors;            Ex=SNeighbors;             ecost=dcOpt.exclusionFactor;
% energy=evaluateFG(double(labeling), Dcost, Lcost, Pw, dcOpt.pairwiseFactor, Ex, ecost, tmpState);
energy=evaluateEnergy(alldpoints, Nhood, ...
    double(labeling), mhs, dcOpt, precomp);
E=energy.value;
bestE=energy.value;

% maximum number of alpha expansion iterations
maxIt=5;
% maxIt=0;
converged=0;
alphaExpIt=0;


% olga_labeling
% pause
% main alpha expansion loop
while ~converged && alphaExpIt<maxIt
    alphaExpIt=alphaExpIt+1;
    thisbest=bestE;
    for expOnAlpha=labelOrder
        %         [~, ~, ~, ~, labeling]=doAlphaExpansion(Dcost, Scost, Lcost, TNeighbors);
        
        
        expLab=expOnAlpha;                      % our current alpha
        
        %prune graph
        mh=[];
        if expLab<outlierLabel, mh=mhs(expLab); end
        [pPoints,pTN, pSN, pDcost, pLabeling,keeppts]= ...
            pruneGraph(alldpoints,labeling,Dcost,TNeighbors,SNeighbors,tmpState2,expLab,outlierLabel,mh,dcOpt,T);
        [npLabels, npPoints]=size(pDcost);

        
        pOlgaLabeling=olga_labeling(keeppts);
        alpInds=find(pLabeling==expLab);         % which nodes are alpha
        notAlphasInd=find(pLabeling~=expLab);    % which are not alpha
        numNotAlphaNodes=numel(notAlphasInd);   % how many not alphas
        
        alphasIndicator=zeros(size(pLabeling));  % put a one on all alphas
        alphasIndicator(alpInds)=1;
        
        newLabeling=pLabeling;
        
        % nothing to do if all nodes are alpha already
        if ~isempty(notAlphasInd)
            
            %%
            % quick check
            
            
            curInds=sub2ind(size(pDcost),double(pLabeling),1:size(pDcost,2));
            % pDcost(curInds);
            
            % set up GCO structure
            [nLabels, nPoints]=size(Dcost);
            
            labeling_=labeling;
            
            % evaluate the enegy before expansion
            Pw=TNeighbors;
            Ex=SNeighbors;
            ecost=dcOpt.exclusionFactor;
            %             energy_=evaluateFG(double(labeling), Dcost, Lcost, Pw, dcOpt.pairwiseFactor, Ex, ecost, tmpState);
            energy_=evaluateEnergy(alldpoints, Nhood, ...
                double(labeling), mhs, dcOpt, precomp);
            E_=energy_.value;
            %             E_
            
            
            %             b=labelcount(labeling,nLabels);
            %             labelUsed=~~b;
            %             lu= labelUsed(1:end-1);
            %
            %             proxcost2=proxcost(lu',:);
            %             proxcost2=proxcost2(:,lu);
            %             Eprox_=sum(proxcost2(:));
            %             E_=E_+Eprox_;
            
            E=E_;
            
            
            %             assert(Eprox_==Eprox2)
            %             pause
            
            
            % energy before expansion
            beforeExpEn=E;
            if beforeExpEn<bestE
                bestE=beforeExpEn;
            end
            
            
            %
            newInds=sub2ind(size(pDcost),expLab*ones(1,size(pDcost,2)),1:size(pDcost,2));
            DcostAlpha=[pDcost(curInds); pDcost(newInds)];
            
            
            
            % which labels are not alpha?
            % unfortunately 'unique' sorts automatically, so unsort here
            [notAlphas, m, n]=unique(pLabeling(notAlphasInd));
            notAlphasUnique=notAlphas;
            [~, sm]=sort(m);
            notAlphas=notAlphas(sm);
            %             notAlphasNoOutlier=setdiff(notAlphas,outlierLabel);
            naoidx=find(notAlphasUnique==outlierLabel);
            notAlphasNoOutlier2=notAlphasUnique;
            if ~isempty(naoidx)
                notAlphasNoOutlier2=notAlphasNoOutlier2([1:naoidx-1 naoidx+1:end]);
            end
            
            notAlphasNoOutlier=notAlphasNoOutlier2;
            %             assert((isempty(notAlphasNoOutlier) && isempty(notAlphasNoOutlier2)) || isequal(notAlphasNoOutlier,notAlphasNoOutlier2));
            
            %             if notAlphas(end)==outlierLabel, notAlphasNoOutlier=notAlphasNoOutlier(1:end-1); end
            %             notAlphas
            %             notAlphasNoOutlier
            
            nAux=numel(notAlphas); % we need nAux auxisiliary nodes
            %             nAux
            
            
            % keep only those that are not alphas
            DcostAlpha=DcostAlpha(:,notAlphasInd);
            DcostAlphaAux=DcostAlpha;
            
            % subtract smaller element
            %             [m mind]=min(DcostAlphaAux);
            %             tmpMat=DcostAlphaAux(sub2ind(size(DcostAlphaAux),mind,1:size(DcostAlphaAux,2)));
            %             DcostAlphaAux=DcostAlphaAux-repmat(tmpMat,2,1);
            
            % add auxiliary unaries
            %             DcostAlphaAux(:,end+1:end+nAux)=[Lcost(notAlphas);zeros(1,nAux)];
            DcostAlphaAux(:,end+1:end+nAux)=[zeros(1,nAux);Lcost(notAlphas)];
            
            Auxmat=zeros(2,0);
            LcostInd=[];
            
            trLabeling=pLabeling(notAlphasInd);
            tmp=numNotAlphaNodes;
            
            % add a pairwise term for each label to the auxiliary variable
            for na=notAlphas
                tmp=tmp+1;
                tmp2=find(trLabeling==na);
                
                % loop implementation
                %                 for nn=1:length(tmp2)
                %                     Auxmat=[Auxmat [tmp2(nn); tmp]];
                %                     LcostInd=[LcostInd tmp-numNotAlphaNodes];
                %                 end
                
                % Matrix implementation
                Auxmat=[Auxmat [tmp2; tmp*ones(1,length(tmp2))]];
                LcostInd=[LcostInd (tmp-numNotAlphaNodes)*ones(1,length(tmp2))];
                
            end
            
            
            
            
            % PAIRWISE SMOOTHNESS (FAST)?
            Dpw=[];
            Dpwi=[]; % indicater whether equal or different
            UnAux=zeros(size(DcostAlphaAux));
            
%             if exist(['constructPairwise.' mexext],'file')
            if dcOpt.mex
                if dcOpt.pairwiseFactor
                    TNeighborsAux=full(pTN+pTN');
                    [UnAux2, Dpwa, Dpwb, Dpwi2]= ...
                        constructPairwise(TNeighborsAux,alphasIndicator,pLabeling,-1);
                    Dpw2=[Dpwa-1;Dpwb-1];
                    Utmp=zeros(size(DcostAlphaAux));
                    Utmp(2,1:numNotAlphaNodes)=-UnAux2*dcOpt.pairwiseFactor;
                    UnAux2=Utmp;
                    
                    % test to matlab implementation
                    %                 assert(isequal(Dpw,Dpw2));  assert(isequal(Dpwi,Dpwi2));  assert(isequal(UnAux,UnAux2));
                    
                    UnAux=UnAux2;
                    Dpw=Dpw2; Dpwi=Dpwi2;
                end
            else
                % PAIRWISE SMOOTHNESS (Matlab, slow)
                if dcOpt.pairwiseFactor
                    TNeighborsAux=pTN+pTN';
                    tmp1=0;
                    for nodeId=notAlphasInd
                        tmp1=tmp1+1;
                        nodesNeighbors=find(TNeighborsAux(nodeId,:));
                        
                        tmp2=0;
                        for nnodeId=nodesNeighbors
                            % neighbor is alpha
                            if ~isempty(find(alpInds==nnodeId, 1))
                                UnAux(2,tmp1)=UnAux(2,tmp1)-dcOpt.pairwiseFactor;
                            else % neighbor is not alpha, keep pairwise connection
                                tmp2=find(notAlphasInd==nnodeId);
                                % always put the one with the lower index first
                                if tmp1<tmp2
                                    
                                    Dpw=[Dpw [tmp1;tmp2]];
                                else
                                    Dpw=[Dpw [tmp2;tmp1]];
                                end
                                
                                
                                if pLabeling(nodeId)==pLabeling(nnodeId) % both are equal, so either leave both or flip both
                                    Dpwi=[Dpwi 1];
                                else
                                    Dpwi=[Dpwi 0];
                                end
                                TNeighborsAux(nodeId,nnodeId)=0;TNeighborsAux(nnodeId,nodeId)=0;
                            end
                        end
                    end
                    Dpw=Dpw-1; % c-code is 0 based
                end
            end
            
            
            % EXCLUSIONS (Supermodular)
            Expw=[];
            Expwi=[]; % indicater whether equal or different
            UnAuxE=zeros(size(DcostAlphaAux));
%             if exist(['constructPairwise.' mexext],'file')
            if dcOpt.mex
                
                % EXCLUSIONS FAST
                if dcOpt.exclusionFactor && expLab~=nLabels
                    SNeighborsAux=full(pSN+pSN');
                    [UnAuxE2, Expwa, Expwb, Expwi2]= ...
                        constructPairwise(SNeighborsAux,alphasIndicator,pLabeling,outlierLabel);
                    Expw2=[Expwa-1;Expwb-1];
                    Utmp=zeros(size(DcostAlphaAux));
                    Utmp(2,1:numNotAlphaNodes)=UnAuxE2*dcOpt.exclusionFactor;
                    UnAuxE2=Utmp;
                    
                    % testing
                    %                 assert(isequal(Expw,Expw2));
                    %                 assert(isequal(Expwi,Expwi));
                    %                 assert(isequal(UnAuxE,UnAuxE2));
                    
                    UnAuxE=UnAuxE2;
                    Expw=Expw2; Expwi=Expwi2;
                end
            else
                % MATLAB IMPLEMENTATION (slow)
                % dont care about exclusion if alpha==outlierLabel
                if dcOpt.exclusionFactor && expLab~=nLabels
                    SNeighborsAux=pSN+pSN';
                    
                    tmp1=0;
                    for nodeId=notAlphasInd
                        tmp1=tmp1+1;
                        nodesNeighbors=find(SNeighborsAux(nodeId,:));
                        
                        tmp2=0;
                        for nnodeId=nodesNeighbors
                            
                            % neighbor is alpha
                            if ~isempty(find(alpInds==nnodeId, 1))
                                UnAuxE(2,tmp1)=UnAuxE(2,tmp1)+dcOpt.exclusionFactor;
                            else % neighbor is not alpha, keep pairwise connection
                                tmp2=find(notAlphasInd==nnodeId);
                                if tmp1<tmp2
                                    Expw=[Expw [tmp1;tmp2]];
                                else
                                    Expw=[Expw [tmp2;tmp1]];
                                end
                                
                                % if both are equal, either leave both or flip both
                                % unless they are both outliers
                                if pLabeling(nodeId)==pLabeling(nnodeId) && pLabeling(nodeId)~=outlierLabel
                                    Expwi=[Expwi 1];
                                else
                                    Expwi=[Expwi 0];
                                end
                                SNeighborsAux(nodeId,nnodeId)=0;
                                SNeighborsAux(nnodeId,nodeId)=0;
                            end
                        end
                    end
                    Expw=Expw-1;
                    
                end
            end
            
            
            
            
            
            
            % now add to unaries
            DcostAlphaAux=DcostAlphaAux+UnAux;
            DcostAlphaAux=DcostAlphaAux+UnAuxE;
            
            %% c code is 0 based
            Auxmatcl=Auxmat-1;
            LcostInd=LcostInd-1;
            %%
            LcostWithProxcost=Lcost;
            
            
            
            % add proximity cost to label cost
            %             pProxCost=proxcost(keeppts',:);pProxCost=pProxCost(:,keeppts);
            if expLab~=outlierLabel
                tmpproxcost=proxcost+proxcost';
                tmpproxcost=tmpproxcost(expLab,:);
                LcostWithProxcost(1:end-1)=LcostWithProxcost(1:end-1)+tmpproxcost;
            end
            
            
            Lcostcl=LcostWithProxcost(notAlphas);
            binLabeling=ones(1,size(DcostAlphaAux,2));
            
            
            if dcOpt.disOpt.initSimple
                %             binLabeling=double(randi(2,1,size(DcostAlphaAux,2))-1); %random
                %             binLabeling=ones(1,size(DcostAlphaAux,2));
                truevars=length(pLabeling);
                GCswitched=find( pOlgaLabeling==expLab);
                binLabeling(GCswitched)=1; % switch some to alpha
                % add a pairwise term for each label to the auxiliary variable
                
                tmp=0;
                for na=notAlphas
                    tmp=tmp+1;
                    tmp2=find(pOlgaLabeling==na); % which ones are alpha in GC solution?
                    tmp3=find(pLabeling==na);
                    
                    %                 na
                    %                 tmp3
                    %                 binLabeling(tmp3)
                    %                 tmp
                    %                 all(pOlgaLabeling==na)
                    %                 truevars+tmp
                    binLabeling(truevars+tmp)=1; % keep unless ...
                    if all(binLabeling(tmp3)) % GC solution replaces all of them
                        binLabeling(truevars+tmp)=0;
                    end
                end
                %             binLabeling(truevars+1:end)=1; % all auxiliaries are on
                %             expLab
                %             size(notAlphas)
                %             size(pLabeling)
                %             size(binLabeling)
                % %             size(trLabeling)
                %             size(pOlgaLabeling)
                %             notAlphas
                %             pLabeling
                %             binLabeling
                % %             trLabeling
                %             pOlgaLabeling
                %             pause
                
            else
                binLabeling=ones(1,size(DcostAlphaAux,2));
%                 binLabeling=zeros(1,size(DcostAlphaAux,2));
%                 binLabeling=double(randi(2,1,size(DcostAlphaAux,2))-1); %random
%                 truevars=length(pLabeling);
%                 binLabeling(truevars+1:end)=1;
            end
            
            binLabeling=ones(1,size(DcostAlphaAux,2));
            pwLcost=proxcost+proxcost';
            pwLcost=pwLcost(:,notAlphasNoOutlier);
            pwLcost=pwLcost(notAlphasNoOutlier',:);
            
            assert(size(pwLcost,1)==length(notAlphasNoOutlier));
            %             pwLcost
            
            % pause
            % save('inferparams.mat', 'dcOpt','Dcost*','Aux*','Lcost*','binLabeling','Dp*','Exp*','notAlphas*','pwLcost','nAux');
            % now that everyhing is done, do inference with openGM
            try
%                 dlmwrite('start_labelin.txt',binLabeling);
%                 pause(.1)

                [Eogm, directLabeling]= ...
                    binaryInference(DcostAlphaAux,Auxmatcl,dcOpt.pairwiseFactor, Lcostcl,LcostInd, ...
                    dcOpt.exclusionFactor,binLabeling,Dpw,Dpwi,Expw,Expwi, notAlphasNoOutlier, pwLcost, nAux, ...
                    dcOpt.disOpt.infParam,0,disAlg);
                
                % % %                 allE=[]; dlab=[];
                % % %                 [allE(1), dlab(1).dl]= ...
                % % %                     binaryInference(DcostAlphaAux,Auxmatcl,dcOpt.pairwiseFactor, Lcostcl,LcostInd, ...
                % % %                     dcOpt.exclusionFactor,binLabeling,Dpw,Dpwi,Expw,Expwi, notAlphasNoOutlier, pwLcost, nAux, ...
                % % %                     [100,0,0],0,1);
                % % %                 [allE(2), dlab(2).dl]= ...
                % % %                     binaryInference(DcostAlphaAux,Auxmatcl,dcOpt.pairwiseFactor, Lcostcl,LcostInd, ...
                % % %                     dcOpt.exclusionFactor,dlab(1).dl,Dpw,Dpwi,Expw,Expwi, notAlphasNoOutlier, pwLcost, nAux, ...
                % % %                     [1,1,1],0,2);
                % % %                 [allE(3), dlab(3).dl]= ...
                % % %                     binaryInference(DcostAlphaAux,Auxmatcl,dcOpt.pairwiseFactor, Lcostcl,LcostInd, ...
                % % %                     dcOpt.exclusionFactor,dlab(1).dl,Dpw,Dpwi,Expw,Expwi, notAlphasNoOutlier, pwLcost, nAux, ...
                % % %                     [100,1e-7,0],0,8);
                % % %                 [allE(4), dlab(4).dl]= ...
                % % %                     binaryInference(DcostAlphaAux,Auxmatcl,dcOpt.pairwiseFactor, Lcostcl,LcostInd, ...
                % % %                     dcOpt.exclusionFactor,dlab(1).dl,Dpw,Dpwi,Expw,Expwi, notAlphasNoOutlier, pwLcost, nAux, ...
                % % %                     [0,2],0,9);
                % %
                % %
                % % %                 pause(1)
                % % %                 [minEn, bestInfAlg]=min(allE);
                % % %                 allE
                % % %                 bestInfAlg
                % % %                 pause
                % % %                 Eogm=minEn;
                % % %                 directLabeling=dlab(bestInfAlg).dl;
                
                
                %                 if ~isempty(GCswitched)
                %                  binLabeling2=zeros(1,size(DcostAlphaAux,2));
                %                  binLabeling2(truevars+1:end)=1; % all auxiliaries are on
                %                 [Eogm2, directLabeling2]= ...
                %                     binaryInference(DcostAlphaAux,Auxmatcl,dcOpt.pairwiseFactor, Lcostcl,LcostInd, ...
                %                     dcOpt.exclusionFactor,binLabeling2,Dpw,Dpwi,Expw,Expwi, notAlphasNoOutlier, pwLcost, nAux, ...
                %                     dcOpt.disOpt.infParam,0,dcOpt.disOpt.alg);
                
                %                 if Eogm~=Eogm2 || ~isequal(directLabeling, directLabeling2)
                %
                %                 GCswitched
                %                 [Eogm Eogm2]
                %                 binLabeling
                %                 binLabeling2
                %                 directLabeling
                %                 directLabeling2
                %                 pause(10)
                %                 end
            catch err
                error('binary inference failed: %s', err.identifier);
            end
            
            if ~all(directLabeling==round(directLabeling))
                directLabeling
            end
            labogm=directLabeling;
            labogm=labogm(1:numNotAlphaNodes);
            %             labogm=labogm+1;
            newLabeling(notAlphasInd(labogm==1))=expLab;
            
            newupLabeling=labeling_;
            newupLabeling(keeppts(notAlphasInd(labogm==1)))=expLab;
            
            % check energy
            %             energy=evaluateFG(double(newupLabeling), Dcost, Lcost, Pw, dcOpt.pairwiseFactor, Ex, ecost, tmpState);
            energy=evaluateEnergy(alldpoints, Nhood, ...
                double(newupLabeling), mhs, dcOpt, precomp);
            E=energy.value;
            
            %             b=labelcount(newupLabeling,nLabels);
            %             labelUsed=~~b;
            %             lu= labelUsed(1:end-1);
            %
            %             proxcost2=proxcost(lu',:);
            %             proxcost2=proxcost2(:,lu);
            %             Eprox=sum(proxcost2(:));
            %             E=E+Eprox;
            
            
            
            
            afterExpEn=E;
            
            % take the new one if it is lower
            minEnLab=newLabeling;
            %             beforeExpEn
            %             afterExpEn
            if beforeExpEn<=afterExpEn
%                                 fprintf('Expansion on %i did not improve energy\n',expLab);
%                                 pause
                newLabeling=pLabeling;
            else
                %                 fprintf('Expansion on %i did improve energy\n',expLab);
                pLabeling=newLabeling;
                labeling=newupLabeling;
                bestE=afterExpEn;
            end
            
            sscnt=sscnt+1;
            if ~mod(sscnt,5) || sscnt==1
                    saveEnergySnapshot(alldpoints, detections, Nhood, newupLabeling, mhs, dcOpt, sceneInfo, 11);
            end
            
        end % if ~isempty(notAlphasInd)
        
    end %for expOnAlpha
    
    
    % we are done if energy did not change after trying to expand on all
    % labels
    if bestE>=thisbest
        converged=1;
    end
end

[~, ~, ~, ~, olga_labeling]=doAlphaExpansion(Dcost, Scost, Lcost, TNeighbors);
% energy_olga=evaluateFG(double(olga_labeling), Dcost, Lcost, Pw, dcOpt.pairwiseFactor, Ex, ecost, tmpState);
energy_olga=evaluateEnergy(alldpoints, Nhood, ...
    double(olga_labeling), mhs, dcOpt, precomp);
if energy.value>energy_olga.value
%     fprintf('graph cuts did better\n');
%         Eogm=energy_olga.value;
%         newupLabeling=double(olga_labeling);
end

logm=newupLabeling;

% clear openGM call to release for writing (compiling)
clear binaryInference