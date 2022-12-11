##### Title: Observed vs expected in regional AFS 
##### Date: 9/29/22
library(ff); library(data.table); library(ggplot2)
library(gridExtra); library(scales); library(MASS)
library(tidyr); library(reshape2);library(ggthemes)
library(gee)

setwd("~/zoellner_research/AFS_project/for_submission/")
source("./scripts/afs_paper_functions.R")

##### 1) Window analysis using observed vs expected 
#1) Compute expected Tajima's D as weighted sum 

#Import 100kb windows with counts of each subtype and observed Tajima's D. Note windows < 50 total sites excluded  
#Script is mst_counts_100kb_windows.R
windows <- fread("./data/all_MST_100kb_counts.txt")
windows$total_sites <- rowSums(windows[,3:98])

#Compute genome wide Tajima's D 
afs_files <- list.files("./data/mst_sfs")

stats_df <- data.frame(matrix(ncol = 5, nrow = 0))
for(file in afs_files){
  print(file)
  mst <- substr(file, 1,7)
  d <- compute_D(paste0("./data/mst_sfs/", file))
  
  afs <- fread(paste0("./data/mst_sfs/", file))
  afs <- afs[-1, ]
  props <- prop.table(afs$V1)
  stats_df <- rbind(stats_df, c(mst, d, props[1], props[2], props[3]))
}

colnames(stats_df) <- c("MST", "D", "singles", "doubles", "triples")
stats_df$D <- as.numeric(stats_df$D)
stats_df$singles <- as.numeric(stats_df$singles); stats_df$doubles <- as.numeric(stats_df$doubles)
stats_df$triples <- as.numeric(stats_df$triples)

windows$D_expect <- 0; windows$singles_expect <- 0
windows$doubles_expect <- 0; windows$triples_expect <- 0

for(i in 1:nrow(windows)){
  print(i)
  windows$D_expect[i] <- weighted.mean(stats_df$D, windows[i, 3:98])
  windows$singles_expect[i] <- weighted.mean(stats_df$singles, windows[i, 3:98])
  windows$doubles_expect[i] <- weighted.mean(stats_df$doubles, windows[i, 3:98])
  windows$triples_expect[i] <- weighted.mean(stats_df$triples, windows[i, 3:98])
}

windows$chr_window <- paste0(windows$chr, ":", windows$window)

windows$orig_diff <- windows$obs_D - windows$D_expect
windows$singles_diff <- windows$singles_prop - windows$singles_expect
windows$doubles_diff <- windows$doubles_prop - windows$doubles_expect
windows$triples_diff <- windows$triples_prop - windows$triples_expect

sd(windows$orig_diff); sd(windows$singles_diff); sd(windows$doubles_diff); sd(windows$triples_diff)

windows$scale_orig_diff <- scale(windows$obs_D - windows$D_expect, center=FALSE)
windows$scale_singles_diff <- scale(windows$singles_prop - windows$singles_expect, center=FALSE)
windows$scale_doubles_diff <- scale(windows$doubles_prop - windows$doubles_expect, center=FALSE)
windows$scale_triples_diff <- scale(windows$triples_prop - windows$triples_expect, center=FALSE)

#Get summary for obs and expected statistics 
mean(windows$obs_D); sd(windows$obs_D)
mean(windows$singles_prop); sd(windows$singles_prop)
mean(windows$doubles_prop); sd(windows$doubles_prop)
mean(windows$triples_prop); sd(windows$triples_prop)

mean(windows$D_expect); sd(windows$D_expect)
mean(windows$singles_expect); sd(windows$singles_expect)
mean(windows$doubles_expect); sd(windows$doubles_expect)
mean(windows$triples_expect); sd(windows$triples_expect)

#2) Plot observed - expected in one plot 
data_long1 <- gather(windows, condition, measurement, scale_orig_diff:scale_triples_diff, factor_key = TRUE)
levels(data_long1$condition) <- c("Tajima's D", "Singleton %", "Doubleton %", "Tripleton %")

#(Figure 6)
ggplot(data_long1, aes(x=measurement, color=(condition))) + geom_density() +
  xlab("Observed - Expected Statistic (Scaled)") + ylab("Density") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color=guide_legend("AFS Statistic")) +
  scale_colour_colorblind() + geom_vline(xintercept = 0, linetype="dotdash", 
                                         color = "blue", size=1)

data_long2 <- gather(windows, condition, measurement, orig_diff:triples_diff, factor_key = TRUE)

#2.5) Add gc content 
gc_content <- read.table("./data/gc_content_100kb.txt", header=FALSE)
colnames(gc_content) <- c("chr","window_start","window_end", "at_pct","gc_pct","num_A","num_c","num_g","num_t","num_N","num_other","seq_length")

gc_content1 <- gc_content[, c(1,2,5)]
colnames(gc_content1) <- c("chr", "window", "gc_pct")
gc_content1$chr <- substr(gc_content1$chr,4,5)
gc_content1$chr_window <- paste0(gc_content1$chr, ":", gc_content1$window)

#3) Regression observed D ~ expected D + RR 
merged_windows <- merge(windows, gc_content1[,c("gc_pct", "chr_window")], by='chr_window')

fit_D <- gee(obs_D ~ D_expect + RR + gc_pct, data = merged_windows, id = chr, corstr = "exchangeable")
D_output <- get_pvals(fit_D)
D_output 

fit_singles <- gee(singles_prop ~ singles_expect + RR + gc_pct, data = merged_windows, id = chr, corstr = "exchangeable")
singles_output <- get_pvals(fit_singles)
singles_output

fit_doubles <- gee(doubles_prop ~ doubles_expect + RR + gc_pct, data = merged_windows, id = chr, corstr = "exchangeable")
doub_output <- get_pvals(fit_doubles)
doub_output

fit_triples <- gee(triples_prop ~ triples_expect + RR + gc_pct, data = merged_windows, id = chr, corstr = "exchangeable")
trip_output <- get_pvals(fit_triples)
trip_output

#4) Investigate negative coefficient for singleton proportions 
has_rr <- subset(merged_windows, is.na(merged_windows$RR) == FALSE)

#(FIGURE S4)
ggplot(data = has_rr, aes(x=singles_expect, y=singles_prop)) + geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("Expected Singleton Proportion") + ylab("Observed Singleton Proportion")
