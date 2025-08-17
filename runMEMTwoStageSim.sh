#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o joblog.$JOB_ID.$TASK_ID #creates a file called joblog.jobidnumber.taskidnumber to write to.
#$ -j y
#$ -l h_rt=20:00:00,h_data=2G #requests 20 hours, 2GB of data (per core)
#$ -pe shared 4 #requests 4 cores
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system
# Notify when
#$ -m ea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)
#$ -t 1-112:1 # 1 to 112, with step size of 1

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load gcc/11.3.0
module load glpk/5.0
module load R/4.5.0-BIO


echo ${SGE_TASK_ID}
echo $1

echo Running MEM_simulation.R for sim_i = ${SGE_TASK_ID} #prints this quote to joblog.jobidnumber
Rscript test/MEM_simulation.R ${SGE_TASK_ID}
