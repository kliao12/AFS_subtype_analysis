#!/bin/bash
#SBATCH --array=1-25
#SBATCH --job-name=sec_fsc
#SBATCH --partition=main
#SBATCH --mem-per-cpu=3Gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<ksliao@umich.edu>
#SBATCH --output=/net/wonderland/home/ksliao/zoellner_research/AFS_project/logs/sec_fsc%a.out

cd /net/wonderland/home/ksliao/zoellner_research/AFS_project/for_submission/data/simulated_neutral

/net/wonderland/home/ksliao/fsc26_linux64/fsc26 -i /net/wonderland/home/ksliao/zoellner_research/AFS_project/for_submission/data/simulated_neutral/mst_par_files/A_C.AAA$SLURM_ARRAY_TASK_ID.par -n100 -q -x -s0 -d -k 2500000 -c8 -r $SLURM_ARRAY_TASK_ID


 
