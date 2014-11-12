#!/bin/tcsh
#
# LSF batch script to run WRF with MPI
#
#BSUB -P P93300075		# project code
#BSUB -W 01:00          	# wall-clock time (hrs:mins)
#BSUB -n 1              	# number of tasks in job
#BSUB -R "span[ptile=1]"	# run 1 MPI task per node
#BSUB -J wrfaer		        # job name
#BSUB -o wrfaer.%J.out          # output file name in which %J is replaced by the job ID
#BSUB -e wrfaer.%J.err	        # error file name in which %J is replaced by the job ID
#BSUB -q caldera		# queue name


# run executable
mpirun.lsf ./wrf.exe


    


