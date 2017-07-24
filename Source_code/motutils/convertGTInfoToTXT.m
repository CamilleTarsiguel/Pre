function allData=convertGTInfoToTXT(gtInfo, txtFile)
% convert Anton's struct format to simple CSV file


gtInfo.Xi=gtInfo.Xi';
gtInfo.Yi=gtInfo.Yi';
gtInfo.W=gtInfo.W';
gtInfo.H=gtInfo.H';

numBoxes = numel(find(gtInfo.Xi(:)));
exGT = find(gtInfo.Xi(:));
[id,fr]=find(gtInfo.Xi);
wd = gtInfo.W(exGT);
ht = gtInfo.H(exGT);
bx = gtInfo.Xi(exGT)-wd/2;
by = gtInfo.Yi(exGT)-ht;

scores=-1*ones(numBoxes,1);

% world coordinates
X=-1*ones(numBoxes,1);
Y=-1*ones(numBoxes,1);
Z=-1*ones(numBoxes,1);

if isfield(gtInfo,'Xgp');
    X=gtInfo.Xgp(exGT);
    Z(:)=0;
end
if isfield(gtInfo,'Ygp');
    Y=gtInfo.Ygp(exGT);
end
if isfield(gtInfo,'Zgp');
    Z=gtInfo.Zgp(exGT);
end

allData = [fr, id, bx, by, wd, ht, scores,X,Y,Z];

if nargin>1
    dlmwrite(txtFile,allData);
end
