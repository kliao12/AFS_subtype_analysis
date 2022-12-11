##### Title: Script to create mutation subtype AFS 
##### Date: 11/13/22
library(data.table)

setwd("~/zoellner_research/AFS_project/for_submission/")

#Import genome wide site info 
data <- fread("./data/gw_sites_info.txt")
colnames(data) 

#Loop over list of 96 subtypes and create SFS using allele counts 
#Export separate files and then aggregate into one dataframe for export
mst_list <- read.table("./data/List_of_MST.txt")

all_afs_df <- as.data.frame(matrix(nrow=7113, ncol = 0))
for(mst in sort(mst_list$V1)){
  print(mst)
  subtype_data <- subset(data, data$V4 == mst)
  ac_factor <- factor(subtype_data$V3, levels = 0:7112)
  afs <- as.data.frame((table(ac_factor)))
  to_export <- afs[,c(2,1)]
  
  #write.table(paste0("./data/mst_sfs/", mst, ".txt", row.names=FALSE, colnames=FALSE, quote=FALSE)
  
  if(mst == 'A_C.AAA'){
    all_afs_df <- cbind(all_afs_df, afs)
    colnames(all_afs_df) <- c("site_num", mst)
  }
  else {
    all_afs_df <- cbind(all_afs_df, afs[,2])
    col_num <- dim(all_afs_df)[2]
    colnames(all_afs_df)[col_num] <- mst
  }
}

#write.table(all_afs_df, "./data/mst_sfs/all_subtypes_afs.txt", row.names = FALSE, col.names = TRUE, quote=FALSE)

