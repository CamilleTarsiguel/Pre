%% sav

%%
inifile=fullfile(confdir,'default.ini');
ini=IniConfig();
ini.ReadFile(inifile);

sec='Parameters';
keys = ini.GetKeys(sec);

% param search dimension
ndim=length(keys);
lval=0; uval=2;

resdone=false(1,maxexper);
for r=1:maxexper
    if exist(sprintf('%s/res_%03d.mat',resdir,r),'file')
        resdone(r)=true;
    end
end
resdone=find(resdone);
ndone=numel(resdone);

allmota=30*ones(1,maxexper);
allmota(resdone)=allmets(resdone,12);

%% determine par a,b,...i
numdir='/home/amilan/research/papers/ongoing/anton-pami/code/paramAnal/';
s1=strfind(setting,'P');s2=strfind(setting,'-');
par=setting(s1+1:s2-1);
motafile=fullfile(numdir,[par '.txt']);
dlmwrite(motafile,allmota);

pname='?';
switch (par)
    case 'a'
        pname='outlier';
    case 'b'
        pname='parsimony';
    case 'c'
        pname='unary';
    case 'd'
        pname='persistence';
    case 'e'
        pname='ang. velocity';
    case 'f'
        pname='lin. velocity';
    case 'g'
        pname='prox. cost';
    case 'h'
        pname='exclusion';
    case 'i'
        pname='pairwise';
end




if ndim==1
    u=1:maxexper;
    X=uval*(u-1)/(maxexper-1);
    
    Z=allmota;
    plot(X,Z);
    xlabel('multiplier');
    ylabel('MOTA [%]');
    
pname='?';
switch (par)
    case 'a'
        pname='outlier';
    case 'b'
        pname='parsimony';
    case 'c'
        pname='unary';
    case 'd'
        pname='persistence';
    case 'e'
        pname='ang. velocity';
    case 'f'
        pname='lin. velocity';
    case 'g'
        pname='prox. cost';
    case 'h'
        pname='exclusion';
    case 'i'
        pname='pairwise';
end    
    legend(pname);
    
elseif ndim==2
    gridX=round(sqrt(maxexper));
    gridY=gridX;
    
    % t=1:maxexper;
    % [u,v]=ind2sub([gridX, gridY],t);
    u=1:gridX;
    v=1:gridY;
    X=uval*(u-1)/(gridX-1);
    Y=uval*(v-1)/(gridX-1);
    
    
    % needs to be transposed
    Z=reshape(allmota,gridX,gridY)';
    surfc(X,Y,Z)
    zlim([-100 100]);
    
    view(2);
    colorbar
    title('MOTA [%]')
    
    pname1='?';
    switch (par(1))
        case 'a', pname1='outlier';
        case 'b', pname1='parsimony';
        case 'c', pname1='unary';
        case 'd', pname1='persistence';
        case 'e', pname1='ang. velocity';
        case 'f', pname1='lin. velocity';
        case 'g', pname1='prox. cost';
        case 'h', pname1='exclusion';
        case 'i', pname1='pairwise';
    end
    pname2='?';
    switch (par(2))
        case 'a', pname2='outlier';
        case 'b', pname2='parsimony';
        case 'c', pname2='unary';
        case 'd', pname2='persistence';
        case 'e', pname2='ang. velocity';
        case 'f', pname2='lin. velocity';
        case 'g', pname2='prox. cost';
        case 'h', pname2='exclusion';
        case 'i', pname2='pairwise';
    end
    
    xlabel(['x ' pname1]);
    ylabel(['x ' pname2]);
end

set(gca,'FontSize',16);