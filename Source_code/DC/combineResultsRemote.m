function mr=combineResultsRemote(jobname,maxexper)

if nargin<2
	maxexper=20;
end

format compact
addpath(genpath('../motutils'));
resdir=sprintf('results/%s',jobname);
disp(resdir)
resfiles=dir(sprintf('%s/res_*.mat',resdir));

dfile=[resdir filesep 'bestres.txt'];
if exist(dfile,'file'), delete(dfile); end
diary(dfile);


if isempty(resfiles)
  warning('nothing there yet!\n')
  return;
end

fprintf('done %d experiments\n',length(resfiles));
oneres=load(sprintf('%s/%s',resdir,resfiles(1).name));
%  allscen=find(mets2d(:,1))';

allmets=[];
for r=1:maxexper
    resfile=sprintf('%s/res_%03d.mat',resdir,r);
    if exist(resfile,'file')
        load(resfile);        
	allmets(r,:)=mean(mets2d(allscen,:),1);
	if howToTrack(allscen(1))
	        allmets(r,:)=mean(mets3d(allscen,:),1);
	end
    end    
end

% best MOTA
allmota=allmets(:,12)';
fprintf('%4d  ',1:maxexper); fprintf('\n');
fprintf('%4.1f  ',allmota); fprintf('\n');
[mm, mr]=max(allmota);


% load best
load(sprintf('%s/res_%03d.mat',resdir,mr));
fprintf('Best: %d\n',mr);

%% write a text file

bestmets=[];
for s=allscen    
    if howToTrack(s)
        tmp=mets3d(s,:);
    else
        tmp=mets2d(s,:);        
    end
    printMetrics(tmp);
    bestmets=[bestmets; tmp];
end

fprintf('--------------------------------------------------------\n');
intmets=4:11;
meanmets=mean(bestmets,1);
meanmets(:,intmets)=round(meanmets(:,intmets));

printMetrics(meanmets);
diary off
