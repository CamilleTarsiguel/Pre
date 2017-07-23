function saveOptimOpts(opt,xmlfile)

% Yangs format always starts with frame 0
% hence subtract the 'correct' first frame

%% Create XML file
docNode = com.mathworks.xml.XMLUtils.createDocument('infOpt');
docRootNode = docNode.getDocumentElement;
% 

%% discrete optimization options
infAlgNames={'TRWS','QPBO','MQPBO','ICM','TRWSi','FastPD','BP','TRBP','LF','LBP'};
algname=char(infAlgNames{opt.disOpt.alg});

disOptNode=docNode.createElement('disOpt');

disOptNode.setAttribute('alg',algname);
if opt.disOpt.alg==1 && opt.disOpt.infParam==3, 
   disOptNode.setAttribute('alg2','LBP');
end

paramsNode=docNode.createElement('infParams');
for p=1:length(opt.disOpt.infParam)
    paramsNode.setAttribute(sprintf('p%d',p),num2str(opt.disOpt.infParam(p)));
end
paramsNode.setAttribute('initSimple',num2str(opt.disOpt.initSimple));
disOptNode.appendChild(paramsNode);

docRootNode.appendChild(disOptNode);

%% continuous optimization options

conOptNode=docNode.createElement('conOpt');

conOptNode.setAttribute('alg',opt.conOpt.alg);

paramsNode=docNode.createElement('infParams');
fn=fieldnames(opt.conOpt);
for f=1:length(fn)
    if ~strcmpi(char(fn{f}),'alg')
        paramsNode.setAttribute(char(fn{f}),num2str(getfield(opt.conOpt,char(fn{f}))));
    end
end
conOptNode.appendChild(paramsNode);

docRootNode.appendChild(conOptNode);

%% save
xmlwrite(xmlfile,docNode);