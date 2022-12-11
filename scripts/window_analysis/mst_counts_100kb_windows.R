##### Title: Script to generate count of each MST, Tajima's D, and prop of singles/doubles/triples 
##### Date: 10/11/22
library(ff); library(data.table); library(ggplot2)
library(gridExtra); library(scales); library(MASS)
library(tidyr)

setwd("~/zoellner_research/AFS_project/for_submission/")

#Import recombination rates
recomb_rates <- read.table("./data/recomb_rate.bed", header = TRUE)

#Import all sites genome wide 
all_data <- fread('./data/gw_sites_info.txt')
colnames(all_data) <- c("chr", "pos", "AC", "MST", "anno")

#Count number of each MST segregating site in window
window_df <- data.frame(matrix(ncol = 103, nrow = 0))

for(i in seq(1,22, by=1)){
    print(i)
    chr_data <- subset(all_data, all_data$chr == i)
    colnames(chr_data) <- c("chr", "pos", "AC", "MST", "anno")

    chr_data_inter <- subset(chr_data, substr(chr_data$anno,1,5) == "Inter")
    chr_data_inter <- subset(chr_data_inter, chr_data_inter$AC > 0)
    chr_recomb_rate <- subset(recomb_rates, CHR == paste0("chr", i))    

    for(lower in seq(min(chr_data_inter$pos), max(chr_data_inter$pos), by = 100000)){
        #print(lower)
        window_data <- subset(chr_data_inter, (lower <= pos & pos <= (lower + 100000)))
        if(nrow(window_data) < 50) {next}
    
        #Get recombination rate for each window 
        chr_window_rr <- subset(chr_recomb_rate, lower <= START  &  END <= (lower+100000))
        recomb_rate <- mean(chr_window_rr$RATE)

        #Compute Tajima's D
        compute_d <- as.data.frame(table(factor(as.numeric(window_data$AC), levels=0:7112)))
        colnames(compute_d) <- c("site_num", "AC_correct")
        compute_d <- subset(compute_d, site_num != 0)

        #Compute pairwise diff
        compute_d$site_num <- as.numeric(as.vector(compute_d$site_num))
        compute_d$n <- 3556*2
        compute_d$pi <- compute_d$AC_correct*compute_d$site_num*(compute_d$n - compute_d$site_num)/ choose(3556*2, 2)
        sum_pi <- sum(compute_d$pi)

        #Compute numerator: Pi - S/a1
        S <- sum(compute_d$AC_correct)
        a1 <- sum(1/seq(1,3556*2)); a2 <- sum(1/seq(1,3556*2)^2)
        n <- 7112

        numer <- sum_pi - S/a1

        #Compute denomenator
        e1 <- (1/a1)*((n+1)/(3*(n-1)) - 1/a1)
        e2 <- 1/(a1^2 + a2)*((2*(n^2 + n + 3)/(9*n*(n-1))) - (n+2)/(n*a1) + (a2/a1^2))

        denom <- sqrt(e1*S + e2*S*(S-1))

        D <- numer/denom

        #Get proportion of singletons/doubletons/tripletons
        singles = prop.table(compute_d$AC_correct)[1]
        doubles = prop.table(compute_d$AC_correct)[2]
        triples = prop.table(compute_d$AC_correct)[3]
        
        #get mst_counts
        mst_counts <- as.data.frame(table(factor(window_data$MST, levels = c("A_C.AAA","A_C.AAC","A_C.AAG","A_C.AAT","A_C.CAA","A_C.CAC","A_C.CAG","A_C.CAT","A_C.GAA","A_C.GAC","A_C.GAG","A_C.GAT","A_C.TAA","A_C.TAC","A_C.TAG","A_C.TAT","A_G.AAA","A_G.AAC","A_G.AAG","A_G.AAT","A_G.CAA","A_G.CAC","A_G.CAG","A_G.CAT","A_G.GAA","A_G.GAC","A_G.GAG","A_G.GAT","A_G.TAA","A_G.TAC","A_G.TAG","A_G.TAT","A_T.AAA","A_T.AAC","A_T.AAG","A_T.AAT","A_T.CAA","A_T.CAC","A_T.CAG","A_T.CAT","A_T.GAA","A_T.GAC","A_T.GAG","A_T.GAT","A_T.TAA","A_T.TAC","A_T.TAG","A_T.TAT","C_A.ACA","C_A.ACC","C_A.ACG","C_A.ACT","C_A.CCA","C_A.CCC","C_A.CCG","C_A.CCT","C_A.GCA","C_A.GCC","C_A.GCG","C_A.GCT","C_A.TCA","C_A.TCC","C_A.TCG","C_A.TCT","C_G.ACA","C_G.ACC","C_G.ACG","C_G.ACT","C_G.CCA","C_G.CCC","C_G.CCG","C_G.CCT","C_G.GCA","C_G.GCC","C_G.GCG","C_G.GCT","C_G.TCA","C_G.TCC","C_G.TCG","C_G.TCT","C_T.ACA","C_T.ACC","C_T.ACG","C_T.ACT","C_T.CCA","C_T.CCC","C_T.CCG","C_T.CCT","C_T.GCA","C_T.GCC","C_T.GCG","C_T.GCT","C_T.TCA","C_T.TCC","C_T.TCG","C_T.TCT"))))
        colnames(mst_counts) <- c("MST", "counts")
        count_wide <- spread(mst_counts, MST, counts)
        window_wide <- cbind(i, lower, count_wide)
        colnames(window_wide)[1:2] <- c("chr", "window")

        if(colnames(window_df)[1] == 'X1'){
            colnames(window_df) <- c(colnames(window_wide), "obs_D", "RR", "singles_prop", "doubles_prop", "triples_prop")  
        }

        window_df[nrow(window_df)+1, ] <- c(window_wide, D,  recomb_rate, singles, doubles, triples)
    }
}

write.table(window_df, paste0("./data/all_MST_100kb_counts.txt"), append = FALSE, sep = " ", dec = ".", col.names=TRUE, row.names=F)