##### Title: Code to compute D-2 statistic for a given frequency spectrum 
##### Date: 12/11/22
##### Author: Kevin Liao (ksliao@umich.edu)

#Function to compute numerator removing singletons and doubletons
compute_newD <- function(mst_afs){
  #Wattersons w/o 1/2 tons
  n <- 3556*2; S <- sum(mst_afs$AC_correct); a1 <- sum(1/seq(1,n-1))
  singles <- mst_afs$AC_correct[1]; doubles <- mst_afs$AC_correct[2]
  theta_w <- (S-singles-doubles)/(a1 - 3/2)
  
  #MPD w/o 1/2 tons 
  mst_afs$pi <- mst_afs$AC_correct*mst_afs$site_num*(n - mst_afs$site_num)
  sum_pi <- sum(mst_afs$pi)
  
  theta_MPD <- (n*(n-1))/((n-2)*(n-3)) * 1/choose(n,2) * (sum_pi - singles*(n-1) - doubles*(n-2)*2)
  newD_num <- theta_MPD - theta_w
  return(newD_num)
}

#Functions to compute denominator using covariance derivations for linear combinations of AFS 
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

afs_cov <- function(c1, c2, sigma, theta, theta_sq){
  an = 0
  for(i in 1:(n-1)){
    an = an + (c1[i]*c2[i]/i)
  }
  
  bn = 0
  for(i in 1:(n-1)){
    for(j in 1:(n-1)){
      bn = bn + c1[i]*sigma[i,j]*c2[j] 
    }
  }
  
  cov = an*theta + bn*theta_sq
  return(cov)
}


####################################################################################################
### Input: dataframe named test_afs with two columns: site number (site_num) and allele count (AC_correct)

# Step 1) Compute numerator 
newD_numer <- compute_newD(test_afs)

# Step 2) Compute Denominator. 
#First compute sigma matrix using results from fu et al. Define haploid sample size 
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

# Define weight vectors for D-2 statistic 
c1 <- vector()
for(i in 1:n){
  c1 <- c(c1, (2*i*(n-i))/((n-2)*(n-3)) )
}
c1[1] = c1[2] = 0

hn <- compute_H(n)
c2 = rep(1/(hn-1.5), n)
c2[1] = c2[2] = 0

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
