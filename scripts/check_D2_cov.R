##### Title: Check compute variance formula 
##### Date: 11/14/22
library(data.table)
library(ggplot2)

setwd("~/zoellner_research/AFS_project/for_submission/")
source("./scripts/afs_paper_functions.R")

#Import precomputed sigma matrix
n=3556*2
sigma <- as.matrix(fread("./output/sigma_mat.txt", header=FALSE))

afs_cov <- function(c1, c2, sigma, theta, theta_sq){
  #get an
  an = 0
  for(i in 1:(n-1)){
    an = an + (c1[i]*c2[i]/i)
  }
  
  #get bn
  bn = 0
  for(i in 1:(n-1)){
    for(j in 1:(n-1)){
      bn = bn + c1[i]*sigma[i,j]*c2[j] 
    }
  }
  
  cov = an*theta + bn*theta_sq
  return(cov)
}

#Define weight vector for Tajima's D
c1 <- vector()
for(i in 1:n){
  c1 <- c(c1, (2*i*(n-i))/(n*(n-1)) )
  #c1 <- c(c1, (2*i*(n-i))/((n-2)*(n-3)) )
}

hn <- compute_H(n)
c2 = rep(1/hn, n)

##### Compute Tajima's D both ways 
compare_D_denoms <- function(subtype){
  temp_AFS <- read.table(paste0("./data/mst_sfs/", subtype, ".txt"))
  temp_AFS <- temp_AFS[-1, ]
  colnames(temp_AFS) <- c("AC_correct", "site_num")
  
  #1) Original way 
  #Compute numerator: Pi - S/a1
  S <- sum(temp_AFS$AC_correct); n <- 7112
  a1 <- sum(1/seq(1,n)); a2 <- sum(1/seq(1,n)^2)

  #Compute denom with original formula 
  e1 <- (1/a1)*((n+1)/(3*(n-1)) - 1/a1)
  e2 <- 1/(a1^2 + a2)*((2*(n^2 + n + 3)/(9*n*(n-1))) - (n+2)/(n*a1) + (a2/a1^2))
  denom <- sqrt(e1*S + e2*S*(S-1))
  
  #Compute denom with weighted approach and sigma matrix 
  Sn = sum(temp_AFS$AC_correct); theta = Sn/hn
  hn_sq <- sum(1/(seq(1:n)^2)); theta_sq = (Sn^2 - Sn)/(hn^2 - hn_sq)
  
  #Get variance of D-2 estimator 
  new_watt_var <- afs_cov(c1, c1, sigma, theta, theta_sq)
  new_mpd_var <- afs_cov(c2, c2, sigma, theta, theta_sq)
  new_mpd_watt_cov <- afs_cov(c1, c2, sigma, theta, theta_sq)
  
  newD_denom <- sqrt(new_mpd_var + new_watt_var - 2*new_mpd_watt_cov)

  return(c(denom, newD_denom))
}

#They're really similar. Could just be small computational differences 
exp1 <- compare_D_denoms("C_T.ACG")
exp2 <- compare_D_denoms("A_C.AAG")




