##### Title: Compute D-2 for each subtype
##### Date: 11/14/22
library(data.table)
setwd("~/zoellner_research/AFS_project/for_submission/")
source("./scripts/afs_paper_functions.R")

args = commandArgs(trailingOnly=TRUE)
subtype <- args[1]
print(subtype)

##### Get variance of D-2 estimator for each subtype
n=3556*2

#Define weight vectors for D-2
c1 <- vector()
for(i in 1:n){
  c1 <- c(c1, (2*i*(n-i))/((n-2)*(n-3)) )
}
c1[1] = c1[2] = 0

hn <- compute_H(n)
c2 = rep(1/(hn-1.5), n)
c2[1] = c2[2] = 0

#Import sigma matrix precomputed 
sigma <- as.matrix(fread("./output/sigma_mat.txt", header=FALSE))

### Run for subtype now ###
test_afs <- read.table(paste0("./data/mst_sfs/", subtype, ".txt"))
test_afs <- test_afs[-1,]
colnames(test_afs) <- c("AC_correct", "site_num")
  
#compute numerator
newD_numer <- compute_newD(test_afs)
  
#compute variance in denominator 
Sn = sum(test_afs$AC_correct)
theta = Sn/hn
  
hn_sq <- sum(1/(seq(1:n)^2))
theta_sq = (Sn^2 - Sn)/(hn^2 - hn_sq)
  
#Get variance of D-2 estimator 
new_watt_var <- afs_cov(c1, c1, sigma, theta, theta_sq)
new_mpd_var <- afs_cov(c2, c2, sigma, theta, theta_sq)
new_mpd_watt_cov <- afs_cov(c1, c2, sigma, theta, theta_sq)
  
newD_denom <- sqrt(new_mpd_var + new_watt_var - 2*new_mpd_watt_cov)
  
#Compute D-2 
newD <- newD_numer/newD_denom

to_output <- as.data.frame(c(subtype, newD))
write.table(to_output, file=paste0("./output/D2_stat/by_subtype/", subtype, ".txt"), row.names = FALSE, col.names = FALSE, quote=FALSE)
