%% installs all necessary dependencies
% Required:
%   - configured mex compiler
%   - internet connection

format compact;

dcdir = pwd;
nerrors=0;

%%%%%%%%%%%%%%%%%%
%%%% MEX Files %%%
%%%%%%%%%%%%%%%%%%
try
    fprintf('Compiling mex...');
    compileMex;
catch err
    fprintf('FAILED: Mex files could not be compiled! %s\n',err.message);
    nerrors=nerrors+1;
end
    



%%%%%%%%%%%%%%%%%%
%%%% MOT Utils %%%
%%%%%%%%%%%%%%%%%%
try
    fprintf('Installing MOT Utils...');
    if ~exist('../motutils','dir')
        cd ..
        !hg clone -v https://bitbucket.org/amilan/motutils
        cd(dcdir);
        fprintf('Success!\n');
    else
        fprintf('IGNORE. Already exist\n');
    end
    
catch err
    fprintf('FAILED: MOT Utils not installed! %s\n',err.message);
    nerrors=nerrors+1;
end

%%%%%%%%%%%%%%%
%%%%%% GCO %%%%
%%%%%%%%%%%%%%%
try
    fprintf('Installing GCO...');
    if ~exist(['external/GCO/matlab/bin/gco_matlab.',mexext],'file')
        cd external
        mkdir('GCO')
        cd GCO
        urlwrite('http://vision.csd.uwo.ca/code/gco-v3.0.zip','gco-v3.0.zip');
        unzip('gco-v3.0.zip');
        cd matlab
        GCO_BuildLib
        
        cd(dcdir);
        fprintf('Success!\n');
    else
        fprintf('IGNORE. Already exist\n');
    end
    
catch err
    fprintf('FAILED: GCO not installed! %s\n',err.message);
    nerrors=nerrors+1;
end

%%%%%%%%%%%%%%%%%%
%%% LIGHTSPEED %%%
%%%%%%%%%%%%%%%%%%
% This one is optional
try
    fprintf('Installing Lightspeed...');
    if ~exist('./external/lightspeed','dir')
        cd external
        mkdir('lightspeed')
        cd lightspeed
        urlwrite('http://ftp.research.microsoft.com/downloads/db1653f0-1308-4b45-b358-d8e1011385a0/lightspeed.zip', ...
            'lightspeed.zip');
        unzip('lightspeed.zip');
        cd lightspeed
        install_lightspeed;
        cd(dcdir);
        
        fprintf('Success!\n');
    else
        fprintf('IGNORE. Already exist\n');
    end
    
catch err
    fprintf('FAILED: Lightspeed not installed! %s\n',err.message);
    nerrors=nerrors+1;
end

%%%%%%%%%%%%%%%%
%%% OpenGM 2 %%%
%%%%%%%%%%%%%%%%
% Only Linux
% Please refer to the [OpenGM website](http://hci.iwr.uni-heidelberg.de/opengm2/) 
% for further instructions. 
% Note that it is possible to compile opengm on Windows, 
% but it is more involved and not officially supported.

if ~ispc
    try
        fprintf('Installing OpenGM...');
        if ~exist('external/opengm','dir')
            cd external
            git clone https://github.com/opengm/opengm.git
            cd opengm
            cd src/external/patches/QPBO
            sh patchQPBO-v1.3.sh
            cd ../TRWS/
            sh patchTRWS-v1.3.sh
            cd ../../../..
            mkdir BUILD
            cd BUILD
            cmake .. -DWITH_QPBO=ON -DWITH_TRWS=ON -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX:PATH=..
            make -j4
            make install
            
            cd(dcdir);
            cd opengm
            compileOGM
            
            cd(dcdir);
            fprintf('Success! Please ...\n');
            fprintf('   ... exit MATLAB,\n');
            fprintf('   ... export library path using ''export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/external/opengm/BUILD/src/external/'' and \n');
            fprintf('   ... start MATLAB again.\n');
        else
            fprintf('IGNORE. Already exist\n');
        end
        
    catch err
        fprintf('FAILED: OpenGM not installed! %s\n',err.message);
        nerrors=nerrors+1;
    end
else
    fprintf('NOTE: OpenGM Installation not supported on Windows\n');
end


fprintf('DC Tracker installed with %d errors\n',nerrors);

cd(dcdir)