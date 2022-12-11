##### Title: Functions for AFS Paper ##### 
##### Date: 9/29/22

#1) compute tajima's D from SFS
compute_D <- function(SFS_path) {
  MST_sfs <-read.table(SFS_path, header = FALSE)
  colnames(MST_sfs) <- c("AC_correct", "site_num")
  MST_sfs <- subset(MST_sfs, site_num != 0)

  #Compute pairwise diff
  MST_sfs$n <- 3556*2
  MST_sfs$pi <- MST_sfs$AC_correct*MST_sfs$site_num*(MST_sfs$n - MST_sfs$site_num)/choose(3556*2, 2)
  sum_pi <- sum(MST_sfs$pi) 

  #Compute numerator: Pi - S/a1
  S <- sum(MST_sfs$AC_correct); n <- 7112
  a1 <- sum(1/seq(1,n)); a2 <- sum(1/seq(1,n)^2)
  numer <- sum_pi - S/a1

  #Compute denomenator
  e1 <- (1/a1)*((n+1)/(3*(n-1)) - 1/a1)
  e2 <- 1/(a1^2 + a2)*((2*(n^2 + n + 3)/(9*n*(n-1))) - (n+2)/(n*a1) + (a2/a1^2))
  denom <- sqrt(e1*S + e2*S*(S-1))
  
  return(numer/denom)
}

#2) Get p values from gee output 
get_pvals <- function(gee_fit){
  output <- as.data.frame(round(summary(gee_fit)$coef, 3))
  output$p_2sided <- 2*pnorm(-abs(output$`Robust z`))
  output$significant <- (output$p_2sided < 0.05/3)
  return(output)
}

#3) Compute new D. Make sure to remove non segregating sites 
compute_newD <- function(mst_afs){
  #Wattersons w/o 1/2 tons
  n <- 3556*2; S <- sum(mst_afs$AC_correct); a1 <- sum(1/seq(1,n-1))
  singles <- mst_afs$AC_correct[1]; doubles <- mst_afs$AC_correct[2]
  theta_w <- (S-singles-doubles)/(a1 - 3/2)
  
  #MPD w/o 1 and 2 tons 
  mst_afs$pi <- mst_afs$AC_correct*mst_afs$site_num*(n - mst_afs$site_num)
  sum_pi <- sum(mst_afs$pi)
  
  theta_MPD <- (n*(n-1))/((n-2)*(n-3)) * 1/choose(n,2) * (sum_pi - singles*(n-1) - doubles*(n-2)*2)
  newD_num <- theta_MPD - theta_w
  return(newD_num)
}

#4) Functions to compute D-2
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

#5) Get T1E rates for null simulations 
get_t1e <- function(subtype, D_df){
  print(subtype)
  data <- D_df[which(D_df[,2] == subtype), ]
  print(head(data))
  tab <- table(abs(data$V1) > 1.96)
  t1E <- tab[2]/(tab[1]+tab[2])
  return(t1E)
}

