##### Title: Script to compute new D on 2500 null simulations for 

library(data.table)
setwd("~/zoellner_research/AFS_project/for_submission/")
source("./scripts/afs_paper_functions.R")

args = commandArgs(trailingOnly=TRUE)
sim_num <- args[1]
subtype <- args[2]

#Lets run this for 2000 neutral simulations across 250 parallel slurm jobs. 10 per job 
chunkit <- split(1:2000, ceiling(seq_along(1:2000)/10))
start <- min(unlist(chunkit[as.numeric(sim_num)]))
end <- max(unlist(chunkit[as.numeric(sim_num)]))

#Import precomputed sigma matrix
n=3556*2
sigma <- as.matrix(fread("./output/sigma_mat.txt", header=FALSE))

#Define weight vectors
c1 <- vector()
for(i in 1:n){
  c1 <- c(c1, (2*i*(n-i))/((n-2)*(n-3)) )
}
c1[1] = c1[2] = 0

hn <- compute_H(n)
c2 = rep(1/(hn-1.5), n)
c2[1] = c2[2] = 0


#Run 5000 neutral sims (this takes a while)
sim_data <- fread(paste0("./data/simulated_neutral/all_", subtype, ".txt"))
#sim_data <- fread(paste0("/net/wonderland/home/ksliao/fsc26_linux64/", subtype, "/", subtype, "_DAFpop0.obs"), skip=2)
#sim_data2 <- sim_data[,c(2:7113)] #get rid of non segregating

sim_newD <- vector()
for(i in start:end){
  print(i)
  temp_AFS <- as.data.frame(cbind(seq(1,7112), t(sim_data[i,])))
  colnames(temp_AFS) <- c("site_num", "AC_correct")
  
  newD_numer <- compute_newD(temp_AFS)
  
  Sn = sum(temp_AFS$AC_correct); theta = Sn/hn
  hn_sq <- sum(1/(seq(1:n)^2)); theta_sq = (Sn^2 - Sn)/(hn^2 - hn_sq)
  
  #Get variance of D-2 estimator 
  new_watt_var <- afs_cov(c1, c1, sigma, theta, theta_sq)
  new_mpd_var <- afs_cov(c2, c2, sigma, theta, theta_sq)
  new_mpd_watt_cov <- afs_cov(c1, c2, sigma, theta, theta_sq)
  
  newD_denom <- sqrt(new_mpd_var + new_watt_var - 2*new_mpd_watt_cov)
  
  #Compute D-2 
  newD <- newD_numer/newD_denom
  sim_newD <- c(sim_newD, newD)
}

sim_D_df <- as.data.frame(sim_newD)

write.table(sim_D_df, paste0("./output/D2_stat/", subtype, "/", sim_num, ".txt"), row.names = FALSE, col.names=FALSE, quote=FALSE) 
