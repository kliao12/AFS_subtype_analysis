##### Title: Script to get all neutral AFS in one file 
##### Date: 11/28/22
library(data.table)

setwd("~/zoellner_research/AFS_project/for_submission/")
source("./scripts/afs_paper_functions.R")

all_sims <- list.files("./data/simulated_neutral/")

AC_df <- as.data.frame(matrix(nrow=0, ncol=7112))
CT_df <- as.data.frame(matrix(nrow=0, ncol=7112))
for(file in all_sims){
  print(file)
  if(substr(file,1,3) == 'A_C'){
    full_path <- paste0("./data/simulated_neutral/", file, "/", file, "_DAFpop0.obs")
    if(file.exists(full_path) == FALSE){next}
    
    data <- fread(full_path, skip=2)
    data <- data[, c(-1, -7114)] #get rid of non segregating and NA column caused by header output
    AC_df <- rbind(AC_df, data)
  } else if(substr(file,1,3) == 'C_T'){
    full_path <- paste0("./data/simulated_neutral/", file, "/", file, "_DAFpop0.obs")
    if(file.exists(full_path) == FALSE){next}
    
    data <- fread(full_path, skip=2)
    data <- data[, c(-1, -7114)] #get rid of non segregating and NA column caused by header output
    CT_df <- rbind(CT_df, data)
  } else{
    next
  }
}

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



write.table(AC_df, "./data/simulated_neutral/all_A_C.AAA.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(CT_df, "./data/simulated_neutral/all_C_T.ACG.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)




### test run 

#Define weight vectors
n=7112
c1 <- vector()
for(i in 1:n){
  c1 <- c(c1, (2*i*(n-i))/((n-2)*(n-3)) )
}
c1[1] = c1[2] = 0

hn <- compute_H(n)
c2 = rep(1/(hn-1.5), n)
c2[1] = c2[2] = 0

#import sigma
sigma <- as.matrix(fread("./output/sigma_mat.txt", header=FALSE))

#import data
sim_data <- fread("./data/simulated_neutral/all_A_C.AAA.txt")
sim_newD <- vector()
for(i in 1:10){
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



