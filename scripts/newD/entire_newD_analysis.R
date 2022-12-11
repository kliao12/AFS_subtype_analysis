##### Title: All D-2 statistic script
##### Date: 11/14/22
setwd("~/zoellner_research/AFS_project/for_submission/")

#Compute sigma matrix to be used for covariance calculations for linear combinations of AFS 
source("./scripts/newD/compute_var.R")

#Compute new D statistic for each subtype AFS using computed sigma matrix
source("./scripts/newD/compute_D2_subtype.R")
source("./scripts/newD/compute_D2_subtype.slurm")

#Simulate neutral AFS to assess new D statistic under null for two subtypes 
source("./scripts/newD/make_fastcoal.R")
"sbatch ./scripts/newD/run_fastcoal_A_C.AAA.slurm"
"sbatch ./scripts/newD/run_fastcoal_C_T.ACG.slurm"

# Look at distribution for neutral simualtions 
#source("./scripts/newD/newD_null_sims.R")
#sbatch ./scripts/newD/newD_null_sims.slurm

null_D_df <- as.data.frame(matrix(nrow=0, ncol=2))
for(subtype in c("A_C.AAA", "C_T.ACG")){
  files <- list.files(paste0("./output/D2_stat/", subtype, "/"), full.names = TRUE)
  
  for(file in files){
    data <- read.table(file, header=FALSE)
    data$subtype = subtype
    null_D_df <- rbind(null_D_df,data)
  }
}

cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data=null_D_df, aes(x=V1, fill=subtype)) + geom_density(alpha=0.5) +
  xlab("D-2 Statistic in Neutral Simulations") + ylab("Density") +
  scale_fill_manual(values=cbp2) + guides(fill=guide_legend(title="Subtype"))

tapply(null_D_df$V1, null_D_df$subtype, summary)    # Summary by group using tapply

#Get T1E for both subtype null simulations 
get_t1e("A_C.AAA", null_D_df)
get_t1e("C_T.ACG", null_D_df)