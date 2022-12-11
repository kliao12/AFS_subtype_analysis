##### Title: Script to run DaDi 10 times for each MST 
##### Date: 11/4/22
setwd("/net/wonderland/home/ksliao/zoellner_research/AFS_project/for_submission/")

#Make sure to run this from conda activate DaDi_3

for(i in 1:10){
  comm1 = paste0('sbatch -J ', i, ' ./scripts/DaDi/run_growth.slurm ', i)
  system(comm1)
  
  comm2 = paste0('sbatch -J ', i, ' ./scripts/DaDi/run_three_epoch.slurm ', i)
  system(comm2)
}

