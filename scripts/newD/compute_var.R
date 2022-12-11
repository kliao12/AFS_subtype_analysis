##### Title: Script to compute variance of D-2 statistic 
##### Date: 11/12/22
library(data.table)
setwd("./zoellner_research/AFS_project/for_submission/")
source("./scripts/afs_paper_functions.R")

#Define functions to compute Hn(i) and Bn(i)
compute_H <- function(n){
  value <- sum(1/seq(1:(n-1)))
  return(value)
}

compute_Bn <- function(i,n){
  first <- (2*n)/((n-i+1)*(n-i))
  second <- compute_H(n+1) - compute_H(i)
  third <- 2/(n-i)
  return(first*second - third)
}

#Compute sigma matrix based off n 
n=3556*2
sigma <- matrix(0, nrow = n-1, ncol = n-1)

#i>j
for(i in 1:(n-1)){
  if(i == 1){next}
  
  for(j in 1:i){
    if(i+j < n){
      sigma[i,j] = (compute_Bn(i+1,n) - compute_Bn(i, n))/2
    }
    
    if(i+j==n){
      first <- (compute_H(n)-compute_H(i))/(n-i)
      second = (compute_H(n)-compute_H(j))/(n-j)
      third = (compute_Bn(i,n) + compute_Bn(j+1,n))/2
      fourth = 1/(i*j) 
      sigma[i,j] = first + second - third - fourth
    }
    
    if(i+j > n){
      first = (compute_Bn(j,n) - compute_Bn(j+1,n))/2
      second = 1/(i*j)
      sigma[i,j] = first - second
    }
  }
}


for(i in 1:(n-1)){
  print(i)

  for(j in 1:i){
    if(i+j < n){
      sigma[j,i] = (compute_Bn(i+1,n) - compute_Bn(i, n))/2
    }
    
    if(i+j==n){
      first <- (compute_H(n)-compute_H(i))/(n-i)
      second = (compute_H(n)-compute_H(j))/(n-j)
      third = (compute_Bn(i,n) + compute_Bn(j+1,n))/2
      fourth = 1/(i*j) 
      sigma[j,i] = first + second - third - fourth
    }
    
    if(i+j > n){
      first = (compute_Bn(j,n) - compute_Bn(j+1,n))/2
      second = 1/(i*j)
      sigma[j,i] = first - second
    }
  }
}

#Set i = j 
for(i in 1:(n-1)){
  if(i < n/2){
    sigma[i,i] = compute_Bn(i+1,n)
  }
  
  if(i == n/2){
    first = 2*(compute_H(n)-compute_H(i))/(n-i)
    second = 1/i^2
    sigma[i,i] = first - second
  }
  
  if(i > n/2){
    sigma[i,i] = compute_Bn(i, n) - 1/i^2
  }
}

#sigma

#Export sigma matrix b/c takes forever to run and can be used across subtypes (sample size n is the same)
#write.table(sigma, file="./output/sigma_mat.txt", row.names=FALSE, col.names=FALSE)


