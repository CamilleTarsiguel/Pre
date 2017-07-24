function allData=convertDetInfoToTXT(detMat, txtFile)
% convert Anton's struct format to simple CSV file


detMat.Xi=detMat.Xi';
detMat.Yi=detMat.Yi';
detMat.W=detMat.W';
detMat.H=detMat.H';
detMat.Sd = detMat.Sd';

numBoxes = numel(find(detMat.Xi(:)));
exDet = find(detMat.Xi(:));
[id,fr]=find(detMat.Xi);
wd = detMat.W(exDet);
ht = detMat.H(exDet);
bx = detMat.Xi(exDet)-wd/2;
by = detMat.Yi(exDet)-ht;

% scores=-1*ones(numBoxes,1);
scores = detMat.Sd(exDet);
id = -1*ones(numBoxes,1);

% world coordinates
X=-1*ones(numBoxes,1);
Y=-1*ones(numBoxes,1);
Z=-1*ones(numBoxes,1);

if isfield(detMat,'Xgp');
    detMat.Xgp=detMat.Xgp';
    X=detMat.Xgp(exDet);
    Z(:)=0;
end
if isfield(detMat,'Ygp');
    detMat.Zgp=detMat.Zgp';
    Y=detMat.Ygp(exDet);
end
if isfield(detMat,'Zgp');
    detMat.Zgp=detMat.Zgp';
    Z=detMat.Zgp(exDet);
end

allData = [fr, id, bx, by, wd, ht, scores,X,Y,Z];

if nargin>1
    dlmwrite(txtFile,allData);
end
