#!/bin/sh
#echo $1

# precomp helper structures 
#matlab -nodesktop -nosplash -r "precomputeForCluster('$PBS_JOBNAME'); quit;";

if [[ $# -eq 2 ]]; then
	echo $2 > config/$1/maxexper.txt
	qsub -t 1-$2 -j oe -l walltime=14:00:00 -l vmem=4GB -o logs -N $1 paramAnal.sh
elif [[ $# -eq 1 ]]; then
        qsub -t 1-20 -j oe -l walltime=14:00:00 -l vmem=4GB -o logs -N $1 paramAnal.sh
else
        echo "Illegal number of parameters"
fi

