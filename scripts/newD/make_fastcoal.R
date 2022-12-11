### Script to simulate neutral SFS for range of theta (Neff, mut rate) ###
library(data.table)
setwd("~/zoellner_research/AFS_project/for_submission/")

source("./scripts/afs_paper_functions.R")

#Import mutation subtype AFS to get mut rates 
merged_data <- read.table("./data/MST_stats_df.txt", header=TRUE)
lambda <- 60/sum(merged_data$nMotifs*merged_data$ERV_rel_rate)
merged_data$abs_mutRate <- merged_data$ERV_rel_rate*lambda

#Get theta estimates for two subtypes. A CpG TpG sites and a random one 
first <- read.table("./data/mst_sfs/A_C.AAA.txt")
first <- first[-1,]
theta_first <- sum(first$V1)/compute_H(3556*2)
mu_first <- subset(merged_data, merged_data$MST == 'A_C.AAA')$abs_mutRate

sec <- read.table("./data/mst_sfs/C_T.ACG.txt")
sec <- sec[-1, ]
theta_sec <- sum(sec$V1)/compute_H(3556*2)
mu_sec <- subset(merged_data, merged_data$MST == 'C_T.ACG')$abs_mutRate

#Need to pass Neff, mu, and length to fastsimcoal
L1 = subset(merged_data, merged_data$MST == 'A_C.AAA')$nMotifs
t1 <- theta_first/L1
Neff_first <- t1/(4*mu_first)

L2 = subset(merged_data, merged_data$MST == 'C_T.ACG')$nMotifs
t2 <- theta_sec/L2
Neff_sec <- theta_sec/(4*mu_sec*L2)

#Create directories for .par files and simulated AFS
system("mkdir ./data/simulated_neutral/mst_par_files")

make_fsc_file <- function(mst, Neff, num_loci, mut_rate, sim_num){
  filename <- paste0("./data/simulated_neutral/mst_par_files/", mst, sim_num, ".par")
  
  cat("//Number of population samples (demes)\n", file = filename, append=FALSE)
  cat("1\n", file = filename, append=TRUE)
  cat("//Population effective sizes (number of genes)\n", file = filename, append=TRUE)
  cat(paste0(Neff, "\n"), file = filename, append=TRUE)
  cat("//Sample sizes\n", file = filename, append=TRUE)
  cat("7112\n", file = filename, append=TRUE)
  cat("//Growth rates  : negative growth implies population expansion\n", file = filename, append=TRUE)
  cat("0\n", file = filename, append=TRUE)
  cat("//Number of migration matrices : 0 implies no migration between demes\n", file = filename, append=TRUE)
  cat("0\n", file = filename, append=TRUE)
  cat("//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix\n", file = filename, append=TRUE)
  cat("0  historical event\n", file = filename, append=TRUE)
  cat("//Number of independent loci [chromosome]\n", file = filename, append=TRUE)
  cat("1 0\n", file = filename, append=TRUE)
  cat("//Per chromosome: Number of linkage blocks\n", file = filename, append=TRUE)
  cat("1\n", file = filename, append=TRUE)
  cat("//per Block: data type, num loci, rec. rate and mut rate + optional parameters\n", file = filename, append=TRUE)
  cat(paste0("DNA ", num_loci, " 0 ", mut_rate),file = filename, append=TRUE)
}

#Run this in parallel across 25 jobs 
for(i in 1:25){
  make_fsc_file("A_C.AAA", round(Neff_first), format(1000000, scientific=FALSE), format(round(mu_first, 10), scientific = FALSE), i)
  make_fsc_file("C_T.ACG", round(Neff_sec), format(1000000, scientific=FALSE), format(round(mu_sec, 10), scientific = FALSE), i)
}
