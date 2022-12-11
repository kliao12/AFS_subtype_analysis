##### Title: Script to generate bootstrapped SFS for each MST via DaDi recs
##### Date: 10/6/22
library(ff); library(data.table); library(ggplot2)
library(gridExtra); library(scales); library(MASS)
library(tidyr)

args <- commandArgs(TRUE)
input_MST <- args[1]
print(input_MST)

setwd("~/zoellner_research/AFS_project/for_submission/")

### Step 1) Construct list of 2Mb windows with their SFS as suggested by https://dadi.readthedocs.io/en/latest/examples/fs_from_data/fs_from_data/
gw_data <- fread("./data/gw_sites_info.txt")
colnames(gw_data) <- c("chr", "pos", "AC", "MST", "anno")

#Running this locally b/c slurm not working 10/6/22
mst_files <- list.files("/net/wonderland/home/ksliao/zoellner_research/AFS_project/for_submission/data/mst_sfs")

for(MST in mst_files){
  print(MST)
  input_MST <- substr(MST, 1, 7)
  
  #make directory
  system(paste0("mkdir ~/zoellner_research/AFS_project/for_submission/data/DaDi/bootstrap_SFS/", input_MST))
  
  gw_MST_data <- subset(gw_data, MST == input_MST)
  
  twoMb_list <-list()
  for(i in seq(1,22, by=1)){
    chr_all_sites <- subset(gw_data, chr == i) #Need this to get same lower/upper for each MST 
    
    chr_data <- subset(gw_MST_data, chr == i) #Get MST for chromosome i
    colnames(chr_data) <- c("chr", "pos", "AC", "MST", "anno")
    
    for(lower in seq(min(chr_all_sites$pos), max(chr_all_sites$pos), by = 2000000)){
      window_id <- paste0("chr", i, "_", lower, ":", lower+2000000)
      window_data <- subset(chr_data, (lower <= pos & pos <= (lower + 2000000)))
  
      window_sfs <- as.data.frame(table(factor(as.numeric(window_data$AC), levels=0:7112)))
      colnames(window_sfs) <- c("site_num", "AC_correct")
      
      twoMb_list[[window_id]] <- window_sfs
    }
  }
  
  
  ### Step 2) Sample with replacement to get bootstrapped datasets 
  num_boot = 100
  for(i in 1:num_boot){
    print(i)
    sampled_list <- sample(twoMb_list, size = 1407, replace = TRUE)
    ac_list <- lapply(sampled_list,"[",2)
    sampled_sfs <- Reduce("+", ac_list) 
    
    write.table(sampled_sfs, paste0("./for_submission/data/DaDi/bootstrap_SFS/", input_MST, "/", input_MST, "_boot_sfs", i, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
  }

}

