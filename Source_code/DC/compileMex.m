% compile mex files
% ...
% 


%% compile utils c
% srcFiles={'labelcount.c','xdotyMex.c','constructPairwise.cpp'};
srcFiles={'labelcount.c','xdotyMex.c', ...
    'EdatBS_mex.c','ElinBS_mex.c','EangBS_mex.c','EexcBS_mex.c','EperBS_mex.c','EbndBS_mex.c'};
srcdir=fullfile('mex','src');
outdir=fullfile('mex','bin');

if ~exist(outdir,'dir'), mkdir(outdir); end

opts='-silent -largeArrayDims CFLAGS="\$CFLAGS -std=c99"'; if ispc, opts=''; end

for k=1:length(srcFiles)
    sprintf('mex %s -outdir %s %s',opts,outdir,fullfile(srcdir,char(srcFiles(k))))
    eval(sprintf('mex %s -outdir %s %s',opts,outdir,fullfile(srcdir,char(srcFiles(k)))));
end

%% compile utils cpp
cppAvailable=~isempty(mex.getCompilerConfigurations('CPP'));
if cppAvailable
    srcFiles={'constructPairwise.cpp'};
    srcdir=fullfile('mex','src');
    outdir=fullfile('mex','bin');

    if ~exist(outdir,'dir'), mkdir(outdir); end

    opts='-silent'; if ispc, opts=''; end

    for k=1:length(srcFiles)
        sprintf('mex %s -outdir %s %s',opts,outdir,fullfile(srcdir,char(srcFiles(k))))
        eval(sprintf('mex %s -outdir %s %s',opts,outdir,fullfile(srcdir,char(srcFiles(k)))));
    end    
end