#!/bin/sh
#
# Sample PBS job script
#
# Copy this script, customize it and submit to PBS with the `qsub''
# commands:
#
# cp pbs-template.sh myjob-pbs.sh
# {emacs|vi} myjob-pbs.sh
# qsub myjob-pbs.sh
#
# If you want your batch job to inherit all your environment variables,
# use the ``V'' switch:
#
# qsub -V myjob-pbs.sh
#
# or uncomment the following line by removing the initial ``###''
#PBS -V 
### Set the job name
#PBS -N Pf6C2_Stnew

### Run in the queue named "default"
### 'batch' is the only queue available on this qchem2 cluster, by default.
#PBS -q batch

### Remove only the three initial "#" characters before #PBS
### in the following lines to enable:
###
### To send email when the job is completed:
### #PBS -m ae
### #PBS -M your@email.address

### Optionally set destination for your program's output
### Specify localhost and an NFS filesystem to prevent file copy errors.
### #PBS -e localhost:$HOME/hpl-test-intel/openmpi.gcc/HPL_JOB.err
### #PBS -o localhost:$HOME/hpl-test-intel/openmpi.gcc/HPL_JOB.log

### Specify the number of cpus for your job.  This example will allocate 16 cores
### using 4 nodes with 4 processes per node.
###
### You MUST specify some number of nodes or TORQUE will fail to load balance.
###
#PBS -l nodes=1:ppn=2

### Tell PBS how much memory you expect to use. Use units of 'b','kb', 'mb' or 'gb'
### #PBS -l mem=256m

### Tell PBS the anticipated run-time for your job, where walltime=HH:MM:SS
#PBS -l walltime=48:00:00

### Switch to the working directory; by default TORQUE launches processes
### from your home directory.
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo The following processors are allocated to this job:
echo `cat $PBS_NODEFILE`

### have wc read from stdin instead of from a file
### ncpus=`wc -l < $PBS_NODEFILE`
### How to run an OpenMPI program
###
### /shared/openmpi/gcc/bin/mpirun -machinefile $PBS_NODEFILE -np $ncpus xhpl

### Or, just run your program
#python3 pdb.main_nongauss.py &
#python3 pdb.main_ngpfine.py &
#python3 pdb.main_vanhove_s.py &
#python3 pdb.main_string.py &
#python3 pdb.main_ht1.py &
#python3 pdb.main_ht2.py &
#python3 pdb.main_ht3.py &
#python3 pdb.main_htlong.py &
python3 pdb.main_St.py &
#python3 pdb.main_Ct.py &
#python3 pdb.main_Ctlong.py &
wait 

### the wait command is a must. the shell needs to wait for all the backgrounded commands.
### without wait, the batch script will quit and kill all the backgrounded jobs.
### see www.olcf.ornl.gov/kb_articles/titan-batch-script-examples/


# These environment variables are available in every batch job
#
# $PBS_ENVIRONMENT set to PBS_BATCH to indicate that the job is a batch job; otherwise,
#                  set to PBS_INTERACTIVE to indicate that the job is a PBS interactive job
# $PBS_JOBID       the job identifier assigned to the job by the batch system
# $PBS_JOBNAME     the job name supplied by the user
# $PBS_NODEFILE    the name of the file that contains the list of nodes assigned to the job
# $PBS_QUEUE       the name of the queue from which the job is executed
# $PBS_O_HOME      value of the HOME variable in the environment in which qsub was executed
# $PBS_O_LANG      value of the LANG variable in the environment in which qsub was executed
# $PBS_O_LOGNAME   value of the LOGNAME variable in the environment in which qsub was executed
# $PBS_O_PATH      value of the PATH variable in the environment in which qsub was executed
# $PBS_O_MAIL      value of the MAIL variable in the environment in which qsub was executed
# $PBS_O_SHELL     value of the SHELL variable in the environment in which qsub was executed
# $PBS_O_TZ        value of the TZ variable in the environment in which qsub was executed
# $PBS_O_HOST      the name of the host upon which the qsub command is running
# $PBS_O_QUEUE     the name of the original queue to which the job was submitted
# $PBS_O_WORKDIR   the absolute path of the current working directory of the qsub command
#
# End of example PBS script

