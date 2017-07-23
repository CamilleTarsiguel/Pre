#!/bin/sh
#echo $1

# precomp helper structures 
#matlab -nodesktop -nosplash -r "precomputeForCluster('$PBS_JOBNAME'); quit;";

qsub -t 1-9 -j oe -l walltime=14:00:00 -l vmem=4GB -o logs -N 1125Pa-1 paramAnal.sh
qsub -t 1-9 -j oe -l walltime=14:00:00 -l vmem=4GB -o logs -N 1125Pb-1 paramAnal.sh
qsub -t 1-9 -j oe -l walltime=14:00:00 -l vmem=4GB -o logs -N 1125Pc-1 paramAnal.sh
qsub -t 1-9 -j oe -l walltime=14:00:00 -l vmem=4GB -o logs -N 1125Pd-1 paramAnal.sh
qsub -t 1-9 -j oe -l walltime=14:00:00 -l vmem=4GB -o logs -N 1125Pe-1 paramAnal.sh
qsub -t 1-9 -j oe -l walltime=14:00:00 -l vmem=4GB -o logs -N 1125Pf-1 paramAnal.sh
qsub -t 1-9 -j oe -l walltime=14:00:00 -l vmem=4GB -o logs -N 1125Pg-1 paramAnal.sh
qsub -t 1-9 -j oe -l walltime=14:00:00 -l vmem=4GB -o logs -N 1125Ph-1 paramAnal.sh
qsub -t 1-9 -j oe -l walltime=14:00:00 -l vmem=4GB -o logs -N 1125Pi-1 paramAnal.sh