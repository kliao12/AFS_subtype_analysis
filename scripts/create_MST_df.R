##### Title: Script to generate dataframe containing all AFS statistics 
##### Date: 11/4/22
library(ggplot2); library(gridExtra); library(scales)
library(MASS); library(data.table); library(gee)
library(ggpubr); library(xlsx)

setwd("/net/wonderland/home/ksliao/zoellner_research/AFS_project/for_submission")
source("./scripts/afs_paper_functions.R")

##### Begin creating dataframe with summary stats for each MST: Tajima D, F*, Ratio of singletons to doubletons, Mutation Rates, etc #####
MST_list <- read.table("./data/List_of_MST.txt")

# Compute Tajimas D using genome wide SFS
taj_d_list = list()

for(i in MST_list$V1){
  if(i == '.'){next}
  
  path_mst = paste0("./data/mst_sfs/", i, ".txt")
  D <- compute_D(path_mst)
  
  taj_d_list[[i]] <- D
  
}

#Create Taj D dataframe from list
taj_d_df <- as.data.frame(cbind(names(taj_d_list), unlist(taj_d_list, use.names=FALSE)))
colnames(taj_d_df) <- c("MST", "D")
taj_d_df$MST <- as.character(taj_d_df$MST)
taj_d_df$D <- as.numeric(as.vector(taj_d_df$D))

# Compute proportion of singletons to doubletons
prop_list = list()
singles_list = list(); doubles_list = list(); triples_list = list(); total_list = list()

for(i in MST_list$V1){
  if(i == '.'){next}
  
  MST_sfs <-read.table(paste0("./data/mst_sfs/", i, ".txt"), header=FALSE)
  colnames(MST_sfs) <- c("AC_correct", "site_num")
  MST_sfs <- subset(MST_sfs, site_num != 0)

  prop_list[[i]] <- MST_sfs[1, "AC_correct"]/MST_sfs[2, "AC_correct"]
  singles_list[[i]] <- prop.table(MST_sfs$AC_correct)[1]
  doubles_list[[i]] <- prop.table(MST_sfs$AC_correct)[2]
  triples_list[[i]] <- prop.table(MST_sfs$AC_correct)[3]
  total_list[[i]] <-   sum(MST_sfs$AC_correct)  
}

props_df <- as.data.frame(cbind(names(singles_list), unlist(singles_list, use.names=FALSE)))
props_df <- cbind(props_df, unlist(doubles_list, use.names=FALSE), unlist(triples_list, use.names=FALSE), unlist(prop_list, use.names=FALSE), unlist(total_list, use.names = FALSE))
colnames(props_df) <- c("MST", "singles", "doubles","triples", "S_D_ratio", "total_variants")

#Import mutation rates
mut_rates <- read.table("./data/ERV_MST_rates.txt", sep=',', header=TRUE)
mut_rates$MST <- ifelse(substr(mut_rates$Type,1,1) == 'A',
                        paste0(substr(mut_rates$Type,1,1), substr(mut_rates$Type,3,4),'.', substr(mut_rates$Motif,1,3)),
                        paste0(substr(mut_rates$Type,2,3), substr(mut_rates$Type,5,5),'.', substr(mut_rates$Motif,1,3))
)

#Create combined dataset
merged_data <- merge(taj_d_df, props_df, by='MST')
merged_data <- merge(merged_data, mut_rates, by='MST')

#write.table(merged_data, "./data/MST_stats_df.txt", col.names = TRUE, row.names = FALSE, quote=FALSE)
#write.xlsx(merged_data, "./data/MST_stats_df.xlsx", row.names = FALSE)

#Subset subtypes to without CpGs
no_CpG <- merged_data[-c(83,87,91,95), ]

CpGs <- c("C_T.ACG","C_T.CCG","C_T.GCG","C_T.TCG")
merged_data$CpG_fill <- ifelse(merged_data$MST %in% CpGs, 2, 1)

merged_data$mt_grouping <- with(merged_data, ifelse(substr(MST,1,3) == 'A_C', 1,0))
merged_data$mt_grouping <- with(merged_data, ifelse(substr(MST,1,3) == 'A_T', 2, merged_data$mt_grouping))
merged_data$mt_grouping <- with(merged_data, ifelse(substr(MST,1,3) == 'C_G', 3, merged_data$mt_grouping))

##### Add on New D-2 estimator 
mst_list <- read.table("./data/List_of_MST.txt")

D2_df <- as.data.frame(matrix(nrow=0, ncol=2))
for(mst in mst_list$V1){
  data <- read.table(paste0("./output/D2_stat/by_subtype/", mst, ".txt"), header=FALSE)
  D2_df <- rbind(D2_df, t(data))
}

colnames(D2_df) <- c("MST", "newD")
D2_df$newD <- as.numeric(D2_df$newD)

merged_data <- merge(merged_data, D2_df, by='MST')


#Write to supplementary table
output_df <- merged_data[, c(1,7,3,4,5,6,2,15,10,11,12)]
write.xlsx(output_df, "./data/MST_stats_df.xlsx", row.names = FALSE)
