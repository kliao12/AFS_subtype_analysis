#!/bin/bash
#SBATCH --array=1-200
#SBATCH --job-name=newD_null
#SBATCH --partition=main
#SBATCH --mem-per-cpu=3Gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<ksliao@umich.edu>
#SBATCH --output=/net/wonderland/home/ksliao/zoellner_research/AFS_project/logs/newD_null%a.out

Rscript /net/wonderland/home/ksliao/zoellner_research/AFS_project/for_submission/scripts/newD_null_sims.R $SLURM_ARRAY_TASK_ID A_C.AAA

Rscript /net/wonderland/home/ksliao/zoellner_research/AFS_project/for_submission/scripts/newD_null_sims.R $SLURM_ARRAY_TASK_ID C_T.ACG
 
