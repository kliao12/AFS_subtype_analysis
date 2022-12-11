##### Title: Script to compute D in 100kb windows and merge on genomic factors 
library(ff); library(data.table); library(ggplot2)
library(gridExtra); library(scales); library(MASS)

#Take as input MST
args <- commandArgs(TRUE)
input_MST <- args[1]

print(input_MST)

#Import recombination rates
recomb_rates <- read.table("/net/wonderland/home/ksliao/zoellner_research/demographic_inf/data/recomb_rate.bed", header = TRUE)

#Import genomewide data
gw_data <- fread("/net/wonderland/home/ksliao/zoellner_research/demographic_inf/data/window_VCFs/genWide_window.txt")
colnames(gw_data) <- c("chr", "pos", "AC", "MST", "anno")

#Rerun analysis on genomewide data by chr
window_df <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(window_df) <- c("chr", "window", "num_weird", "num_total", "prop_weird", "taj_D", "F", "window_RR")

for(i in seq(1,22, by=1)){
    print(i)

    chr_data <- subset(gw_data, chr == i)
    #chr_data <- read.table(paste0("/net/wonderland/home/ksliao/zoellner_research/demographic_inf/data/window_VCFs/",i,"_window.vcf"), header=FALSE)
    colnames(chr_data) <- c("chr", "pos", "AC", "MST", "anno")

    chr_data_inter <- subset(chr_data, substr(chr_data$anno,1,5) == "Inter")
    chr_recomb_rate <- subset(recomb_rates, CHR == paste0("chr", i))
        
    for(lower in seq(min(chr_data_inter$pos), max(chr_data_inter$pos), by = 25000)){
        print(lower)
        window_data <- subset(chr_data_inter, (lower <= pos & pos <= (lower + 25000)))
        if(nrow(window_data) < 50) {next}
  
        chr_window_rr <- subset(chr_recomb_rate, lower <= START  &  END <= (lower+25000))
        recomb_rate <- mean(chr_window_rr$RATE)
  
        mst_counts <- as.data.frame(table(window_data$MST))
        colnames(mst_counts) <- c("MST", "counts")
        weird_counts <- subset(mst_counts, MST == input_MST)
        test_prop <- sum(weird_counts$counts)/sum(mst_counts$counts)
  
        compute_d <- as.data.frame(table(factor(as.numeric(window_data$AC), levels=0:7112)))
        colnames(compute_d) <- c("site_num", "AC_correct")
        compute_d <- subset(compute_d, site_num != 0)  

        #Compute pairwise diff
        compute_d$site_num <- as.numeric(as.vector(compute_d$site_num))
        compute_d$n <- 3556*2
        compute_d$pi <- compute_d$AC_correct*compute_d$site_num*(compute_d$n - compute_d$site_num)/
            choose(3556*2, 2)
        sum_pi <- sum(compute_d$pi) 
              
        #Compute numerator: Pi - S/a1
        S <- sum(compute_d$AC_correct)
        a1 <- sum(1/seq(1,3556*2))
        a2 <- sum(1/seq(1,3556*2)^2)
        n <- 7112
          
        numer <- sum_pi - S/a1
  
        #Compute denomenator
        e1 <- (1/a1)*((n+1)/(3*(n-1)) - 1/a1)
        e2 <- 1/(a1^2 + a2)*((2*(n^2 + n + 3)/(9*n*(n-1))) - (n+2)/(n*a1) + (a2/a1^2))
          
        denom <- sqrt(e1*S + e2*S*(S-1))
          
        D <- numer/denom
          
        #Compute F* while I'm at it
        compute_f <- as.data.frame(table(factor(as.numeric(window_data$AC), levels=0:7112)))
        colnames(compute_f) <- c("site_num", "AC_correct")
        compute_f <- subset(compute_f, site_num != 0)
        sapply(compute_f, typeof)
          
        compute_f$site_num <- as.numeric(as.vector(compute_f$site_num))
        compute_f$AC_correct <- as.numeric(compute_f$AC_correct)

        compute_f$n <- 7112
        compute_f$pi <- compute_f$AC_correct*compute_f$site_num*(compute_f$n - compute_f$site_num)/
          choose(7112, 2)
        sum_pi <- sum(compute_f$pi) 
          
        eps_1 <- compute_f[1, "AC_correct"]
        eps_n_1 <- compute_f[nrow(compute_f), "AC_correct"]
        delta_1_n_1 <- 0
        n <- 7112
        eta_1 <- (eps_1 + eps_n_1) / (1 + delta_1_n_1)
          
        f_numer <- sum_pi - ((n-1)/n)*eta_1
          
        #Compute denomenator: Use formula from Fu and Li's paper
        eta <- sum(compute_f$AC_correct)
        a_n <- sum(1/seq(1,n-1))
        a_n_1 <- sum(1/seq(1,n))
        b_n <- sum(1/seq(1,n-1)^2)
        c_n <- 2*(n*a_n - 2*(n-1))/((n-1)*(n-2))
        d_n <- c_n + (n-2)/(n-1)^2 + (2/(n-1))*(1.5 - (2*a_n_1 - 3)/(n-2) - 1/n)
          
        v_f <- (d_n + 2*(n^2+n+3)/(9*n*(n-1)) - 2*(4*b_n - 6 + 8/n)/(n-1)) / (a_n^2+b_n)
        u_f <- (n/(n+1) + (n+1)/(3*(n-1)) - 4/(n*(n-1)) + (2*(n+1)/(n-1)^2)*(a_n_1 - 2*n/(n+1)))/a_n - v_f
          
        var_f <- u_f*eta + v_f*eta^2
          
        f_denom <- sqrt(var_f)
          
        f <- f_numer / f_denom
          
        window_df[nrow(window_df)+1, ] <- c(i, lower, sum(weird_counts$counts), sum(mst_counts$counts), test_prop, D, f, recomb_rate)
    }

}

window_df1 <- subset(window_df, num_weird > 0)

write.table(window_df1, paste0("/net/wonderland/home/ksliao/zoellner_research/demographic_inf/output/tajD_compare_25000/original/tajD_", input_MST, ".txt"), append = FALSE, sep = " ", dec = ".", col.names=TRUE, row.names=FALSE, quote=FALSE)



