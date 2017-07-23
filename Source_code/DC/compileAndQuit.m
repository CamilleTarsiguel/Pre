function compileAndQuit(fileToCompile)
% compile with MCC

fprintf('compiling %s...\n',fileToCompile);
[pt fl ex]=fileparts(fileToCompile);

compstr=sprintf('mcc -a /gris/gris-f/home/aandriye/software/lightspeed/xrepmat.m -I ~/visinf/projects/ongoing/dctracking/PAMI/opengm  -I ~/visinf/projects/ongoing/dctracking/PAMI/fminlbfgs -I ~/software/gco-v3.0/matlab -I ~/software/gco-v3.0/matlab/bin -I ~/software/lightspeed -I ../../contracking/utils -I ../../contracking/utils/splinefit -I ../../contracking/utils/camera  -I ./mex/bin -I ../../contracking/utils/mex/bin -R -singleCompThread -R -nodisplay -C %s%s -m %s;',fl,ex,fl);
compstr

while true
    try
	fprintf('trying...\n');
	eval(compstr);
	quit;
    catch err
        fprintf('not this time\n');
        pause(15);
    end
end