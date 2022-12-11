##### Title: Look at subtype abundances by Tajima's D Quantiles 
##### Date: 11/14/22
library(gtools); library(tidyr); library(dplyr);
library(ggplot2); library(ggpattern)
setwd("~/zoellner_research/AFS_project/for_submission/")

#Import 100kb windows with tajD computed from any of the 96 subtypes
subtype = 'C_T.GCC' #This is random, doesn't matter
mst_data = read.table("./data/100kb_windows/tajD_C_T.GCC.txt", header=TRUE)
mst_data$tajD_group <- quantcut(mst_data$taj_D, q=20)
top_windows_df <- mst_data[,c(1,2,9)]

mst_list <- read.table("./data/List_of_MST.txt")

#For each of 96 subtypes, add classifier for each window if in top 10% of subtype abundance 
for(i in mst_list$V1){
  print(i)
  if(i == '.'){next}
  
  mst_data <- read.table(paste0("./data/100kb_windows/tajD_", i, ".txt"), header=TRUE)
  mst_data$tajD_sig <- ifelse(mst_data$taj_D < quantile(mst_data$taj_D, probs = 0.05), 1, 0)
  mst_data$prop_weird_sig <- ifelse(mst_data$prop_weird > quantile(mst_data$prop_weird, probs = 0.90), 1, 0)
  keep <- mst_data[, c(1,2,10)]
  colnames(keep) <- c('chr', 'window', i)
  
  top_windows_df <- merge(top_windows_df, keep, by=c('chr','window'), all.x=TRUE)
}

top_windows_df[is.na(top_windows_df)] <- 0

#Check whether mutation rate or gBGC causes subtype abundance to vary across tajima's D quantiles 
merged_data <- read.table("./data/MST_stats_df.txt", header=TRUE)
median_mut_rate <- median(merged_data$ERV_rel_rate)

no_CpG <- merged_data[-c(83,87,91,95), ]

SW_types <- c("C_T", "C_A"); WS_types <- c("A_C", "A_G"); neutral_types <- c("C_G", "A_T")
CpG_TpG <- c("C_T.ACG", "C_T.CCG", "C_T.GCG", "C_T.TCG")

SW_low_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% SW_types & no_CpG$ERV_rel_rate <= median_mut_rate)$MST
SW_high_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% SW_types & no_CpG$ERV_rel_rate > median_mut_rate)$MST

WS_low_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% WS_types & no_CpG$ERV_rel_rate <= median_mut_rate)$MST
WS_high_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% WS_types & no_CpG$ERV_rel_rate > median_mut_rate)$MST

neutral_low_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% neutral_types & no_CpG$ERV_rel_rate <= median_mut_rate)$MST
neutral_high_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% neutral_types & no_CpG$ERV_rel_rate > median_mut_rate)$MST

top_windows_df$CpG_TpG <- rowSums(top_windows_df[,CpG_TpG])
top_windows_df$SW_low_mut <- rowSums(top_windows_df[,SW_low_mut])
top_windows_df$SW_high_mut <- rowSums(top_windows_df[,SW_high_mut])
top_windows_df$WS_low_mut <- rowSums(top_windows_df[,WS_low_mut])
top_windows_df$WS_high_mut <- rowSums(top_windows_df[,WS_high_mut])
top_windows_df$neutral_low_mut <- rowSums(top_windows_df[,neutral_low_mut])
top_windows_df$neutral_high_mut <- rowSums(top_windows_df[,neutral_high_mut])

CpG_TpG_quantiles <- as.data.frame(tapply(top_windows_df$CpG_TpG, top_windows_df$tajD_group, mean))
SW_low_mut_quantiles <- as.data.frame(tapply(top_windows_df$SW_low_mut, top_windows_df$tajD_group, mean))
SW_high_mut_quantiles <- as.data.frame(tapply(top_windows_df$SW_high_mut, top_windows_df$tajD_group, mean))
WS_low_mut_quantiles <- as.data.frame(tapply(top_windows_df$WS_low_mut, top_windows_df$tajD_group, mean))
WS_high_mut_quantiles <- as.data.frame(tapply(top_windows_df$WS_high_mut, top_windows_df$tajD_group, mean))
neutral_low_mut_quantiles <- as.data.frame(tapply(top_windows_df$neutral_low_mut, top_windows_df$tajD_group, mean))
neutral_high_mut_quantiles <- as.data.frame(tapply(top_windows_df$neutral_high_mut, top_windows_df$tajD_group, mean))

all_categories <- cbind(CpG_TpG_quantiles, SW_low_mut_quantiles, SW_high_mut_quantiles,
                        WS_low_mut_quantiles, WS_high_mut_quantiles, neutral_low_mut_quantiles, neutral_high_mut_quantiles)
all_categories$D_quantile <- as.factor(rownames(all_categories))
all_categories$D_quantile_rank <- rank(all_categories$D_quantile)
colnames(all_categories) <- c("CpG_TpG", "SW_low", "SW_high", "WS_low", "WS_high", "Neutral_low", "Neutral_high", "D_quantiles", "Quantile_rank")

categories_long <- gather(all_categories, group, num_abundant, CpG_TpG:Neutral_high, factor_key=TRUE)
categories_long$gBGC <- sub("\\_.*", "", categories_long$group)
categories_long$mut_rate <- sub(".*\\_", "", categories_long$group)
categories_long$gBGC <- ifelse(categories_long$gBGC == 'CpG', 'CpG TpG', categories_long$gBGC)
categories_long$mut_rate <- ifelse(categories_long$mut_rate == 'TpG', 'CpG TpG', categories_long$mut_rate)

#(FIGURE 5)
ggplot(data=categories_long, aes(x=D_quantiles, y=num_abundant, fill=gBGC, pattern=mut_rate)) +
  geom_bar_pattern(position='stack', stat='identity',
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025) + 
  scale_pattern_manual("Mutation Rate", values = c(high = "stripe", low = "none", `CpG TpG`='circle'),
                       labels=c("CpG TpG", "High", "Low")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  xlab("D Quantiles") + ylab("Average # Abundant Subtypes") +
  scale_x_discrete(limits = rev) + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#(FIGURE S3)
all_categories$total_abun <- rowSums(all_categories[,1:7])
all_categories_prop <- all_categories[,1:7]/rowSums(all_categories[,1:7])
all_categories_prop$tajD_group <- rownames(all_categories)
all_categories_prop$D_quantile_rank <- rank(desc(all_categories_prop$tajD_group))

cat_only <- all_categories_prop[,1:7]
(cat_only[nrow(cat_only),] - cat_only[1,]) / cat_only[1,]

prop_long <- gather(all_categories_prop, cat_group, prop_abundant, CpG_TpG:Neutral_high, factor_key=TRUE)

ggplot(data=prop_long, aes(x=reorder(tajD_group, desc(tajD_group)), y=prop_abundant, color=cat_group)) + 
  geom_line(aes(group=cat_group)) +
  geom_point() +
  xlab("D Quantiles") + ylab("Proportion of Abundant Subtypes") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(name = "Category Group")

