##### Title: Script for all DaDi Analysis 
##### Date: 11/3/22
##### Input: SFS for each MST 
library(data.table); library(openxlsx); library(dplyr)
library(ggplot2)

setwd("/net/wonderland/home/ksliao/zoellner_research/AFS_project/for_submission")

### 1) Run DaDi for growth and three epoch model using python script run on slurm cluster 
source("./scripts/DaDi/bootstrap_SFS.R") #this takes a while. needed to get bootstrap SFS for dadi uncertainty analysis
source("./scripts/DaDi/run_dadi_both_models_10times.R") #this runs dadi using run_growth.py/slurm and run_three_epoch.py/slurm

### 2) Combine inferred output into one file for appendix table
mst_list <- read.table("./data/List_of_MST.txt")

growth_df <- as.data.frame(matrix(ncol=8, nrow=0))
three_epoch_df <- as.data.frame(matrix(ncol=12, nrow=0))

for(mst in mst_list$V1){
  for(i in 1:10){
    if(file.exists(paste0("./output/DaDi/", mst, "_growth", i, ".txt")) == FALSE){next}
    if(file.exists(paste0("./output/DaDi/", mst, "_three_epoch", i, ".txt")) == FALSE){next}
    
    g_data <- fread(paste0("./output/DaDi/", mst, "_growth", i, ".txt"))
    t_data <- fread(paste0("./output/DaDi/", mst, "_three_epoch", i, ".txt"))
    
    growth_df <- rbind(growth_df, g_data)
    three_epoch_df <- rbind(three_epoch_df, t_data)
  }
}

growth_df <- growth_df[,-1]
three_epoch_df <- three_epoch_df[,-1]

#Find best fit model for each subtype 
best_growth <- growth_df %>% 
  group_by(MST) %>% 
  slice_max(order_by = ll_model)

best_three_epoch <- three_epoch_df %>% 
  group_by(MST) %>% 
  slice_max(order_by = log_lik)

#Export to txt
#write.table(best_growth, "./output/DaDi/best_MST_growth.txt", row.names = FALSE, col.names = TRUE, quote=FALSE)
#write.table(best_three_epoch, "./output/DaDi/best_MST_three_epoch.txt", row.names = FALSE, col.names = TRUE, quote=FALSE)

### 3) Analyze output and make plots
mut_rates <- read.table("./data/ERV_MST_rates.txt", sep=',', header=TRUE)
mut_rates$MST <- ifelse(substr(mut_rates$Type,1,1) == 'A',
                        paste0(substr(mut_rates$Type,1,1), substr(mut_rates$Type,3,4),'.', substr(mut_rates$Motif,1,3)),
                        paste0(substr(mut_rates$Type,2,3), substr(mut_rates$Type,5,5),'.', substr(mut_rates$Motif,1,3))
)

growth <-read.table("./output/DaDi/best_MST_growth.txt", header=TRUE)
growth1 <- merge(growth, mut_rates, by="MST")

##### Compute absolute mutation rate b/c previous paper used relative rates
lambda <- 60/sum(growth1$nMotifs*growth1$ERV_rel_rate)
growth1$abs_mutRate <- growth1$ERV_rel_rate*lambda

growth1$N_ref <- growth1$theta/(4*growth1$nMotifs*growth1$abs_mutRate)
growth1$T_gen <- 2*growth1$N_ref * growth1$T

#Export for appendix 
#write.xlsx(growth1, "./output/DaDi/best_MST_growth_parameters.xlsx")

#Make graph of singleton proportion vs growth ratio 
mst_df <- read.table("./data/MST_stats_df.txt", header=TRUE)
growth2 <- merge(growth1, mst_df[,c("MST","singles")], by='MST')

growth2$CpG <- as.factor(ifelse(substr(growth2$MST, 1,3) == 'C_T' & substr(growth2$MST,6,7) == 'CG', 2, 1))

CpGs <- c("C_T.ACG","C_T.CCG","C_T.GCG","C_T.TCG")

#(Figure 3)
ggplot(data=growth2, aes(x=singles, y=nu, color = CpG, shape=CpG)) + geom_point(size = 2.5) + 
  xlab("Proportion of Singletons") +   ylab("Current Size / Ancestral Size") +
  theme(axis.text.x=element_text(angle=90)) +
  theme(legend.position = c(0.15, 0.8), legend.title = element_blank(), axis.title.y =element_text(size=12), axis.title.x =element_text(size=12)) +
  scale_color_manual(labels = c(" Non X[C->T]G", " X[C->T]G"),  values=c("royalblue3", "firebrick4")) +
  scale_shape_manual(labels = c(" Non X[C->T]G", " X[C->T]G"),  values=c(19,17))

cor.test(growth2$singles, growth2$nu)

##### Repeat for three_epoch 
three_epoch <- fread("./output/DaDi/best_MST_three_epoch.txt")
three_epoch$time <- three_epoch$TB + three_epoch$TF

MST_stats_df <- read.table("./data/MST_stats_df.txt", header=TRUE)

three_epoch1 <- merge(three_epoch, MST_stats_df, by='MST')

lambda <- 60/sum(three_epoch1$nMotifs*three_epoch1$ERV_rel_rate)
three_epoch1$abs_mutRate <- three_epoch1$ERV_rel_rate*lambda
three_epoch1$N_ref <- three_epoch1$theta/(4*three_epoch1$nMotifs*three_epoch1$abs_mutRate)
three_epoch1$T_gen <- 2*three_epoch1$N_ref * three_epoch1$time

#export for paper
#write.xlsx(three_epoch1, "./output/DaDi/best_MST_three_epoch_parameters.xlsx")

ggplot(data=three_epoch1, aes(x=D, y=nuB)) + geom_point() +
  ggtitle("Population Contraction During Bottleneck by Tajima's D for 96 MSTs") +
  xlab("Tajima's D") + ylab("Population Contraction During Bottleneck") +
  geom_smooth(method='lm') + theme(plot.title = element_text(hjust = 0.5, size = 13))

#Look at correlation
cor.test(three_epoch1$D, three_epoch1$nuB)

#Look at how many MST suggest no bottleneck vs a deep bottleneck
no_bottle <- subset(three_epoch1, three_epoch1$nuB > 0.85)
alot_bottle <- subset(three_epoch1, three_epoch1$nuB < 0.15)

