#!/bin/sh
/usr/local/matlab-2012a/bin/matlab -r "compileAndQuit('$1')"
cp myrun.sh run_randomParamSearchCluster.sh
