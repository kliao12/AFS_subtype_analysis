###### Title: Script to run all analysis in The Effect of Mutation Subtypes on Population Genetics Inference 
###### Date: 11/11/22
###### Author: Kevin Liao 
setwd("~/zoellner_research/AFS_project/for_submission/")
source("./scripts/afs_paper_functions.R")
source("./scripts/libraries_needed.R")

### 0) Create AFS for each of 96 mutation subtypes. Outputs separate txt file for each and a single with all 96 subtypes.
#source("./scripts/create_subtype_AFS.R")

### 0.5) Analysis for new D-2 statistic. Note some of these steps take a while running from scratch. 
### Precomputed D-2 values can be imported in step 1)
#source("./scripts/newD/entire_newD_analysis.R")

### 1) Create dataframe of genome wide AFS statistics for each of 96 subtypes. Creates df called merged_data
source("./scripts/create_MST_df.R")

#(FIGURE S1): Create plot for Tajimas D by MST
taj_d <- ggplot(data=merged_data, aes(x=MST, y=D, fill = as.factor(mt_grouping)) ) + 
  geom_bar(stat="identity", show.legend = FALSE) + 
  xlab("Mutation Subtype") + ylab("D") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=90)) +
  scale_y_continuous(limits=c(-2.3, -1.4),oob = rescale_none) +
  scale_fill_manual(values = c("3"="cornflowerblue", "2" = "cornflowerblue",
                               "1" = "cornflowerblue", "0" = "grey70"))
taj_d

### 2) AFS heterogeneity across subtypes and signals of recurrent mutations and gene conversion shaping the AFS 

#1) Signals of mutation rate/recurrent mutations on AFS 
#Correlation between mutation rate and Singleton to Doubleton Ratio
cor.test(no_CpG$ERV_rel_rate, no_CpG$S_D_ratio); cor.test(merged_data$ERV_rel_rate, merged_data$S_D_ratio)

merged_data$MT <- substr(merged_data$MST, 1, 3)
merged_data$MT_arrow <- paste0(substr(merged_data$MT,1,1), "->", substr(merged_data$MT,3,3))

#(FIGURE 3): Plot mutation rate vs SD ratio overall and broken down by 6 mutation types 
a1 <- ggplot(merged_data, aes(x=ERV_rel_rate, y=S_D_ratio)) +
  geom_point(aes(color=as.factor(CpG_fill))) + 
  xlab("Mutation Rate") + ylab("Singletons / Doubletons") + geom_smooth(method = "lm", se = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits=c(3.5,8),oob = rescale_none) +
  scale_color_manual(values = c("2" = "orangered2", "1" = "black"), name="Mutation\nSubtype", breaks = c(1,2), labels =c("Other","CpG TpG")) 
a1

equal_breaks <- function(n = 4, s = 0.05, ...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq(min(x)+d, max(x)-d, length=n)
  }
}

all_six <- ggplot(data=merged_data, aes(x=ERV_rel_rate, y=S_D_ratio)) + geom_point(aes(color=as.factor(CpG_fill))) + geom_smooth(method = "lm", se=FALSE) + 
  facet_wrap(~MT_arrow, scales = "free_x") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Mutation Rate") + ylab("Singletons / Doubletons") + 
  scale_color_manual(values = c("2" = "orangered2", "1" = "black"),
                     name="Mutation\nSubtype",
                     breaks = c(1,2),
                     labels =c("Other","CpG TpG")) +
  theme(panel.spacing = unit(2, "lines")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.0001), breaks=equal_breaks(n=4))
all_six
ggarrange(a1, all_six, ncol=1,labels=c("a)", "b)"))

#2) Signals of gBGC on AFS. Split into WS, SW, and indifferent mutation types. Compare D and D-2 across subtypes  
#Remove CpGs as done in Lachance et al (AJHG, 2014)
remove_CpG <- subset(merged_data, substr(merged_data$MST,6,7) != 'CG')

WS <- subset(remove_CpG, substr(MST,1,3) %in% c("A_C", "A_G"))
non_WS <- subset(remove_CpG, !(substr(MST,1,3) %in% c("A_C", "A_G")))
t.test(WS$D, non_WS$D)

SW <- subset(remove_CpG, substr(MST,1,3) %in% c("C_A", "C_T"))
non_SW <- subset(remove_CpG, !(substr(MST,1,3) %in% c("C_A", "C_T")))
t.test(SW$D, non_SW$D)

neutral <- subset(remove_CpG, substr(MST,1,3) %in% c("C_G","A_T"))
non_neutral <- subset(remove_CpG, !(substr(MST,1,3) %in% c("C_G","A_T")))
t.test(neutral$D, non_neutral$D)

#Repeat for new D-2 statistic  
t.test(WS$newD, non_WS$newD)
t.test(SW$newD, non_SW$newD)
t.test(neutral$newD, non_neutral$newD)

### 3) Effect of genome wide AFS heterogeneity across subtypes on demographic inference using DaDi
source("./scripts/DaDi/all_DaDi_analysis.R")

### 4) Effect of subtype composition and local genomic factors on the local AFS 
#source("./scripts/window_analysis/mst_counts_100kb_windows.R") #Run to get df of 100kb windows w/ counts of 96 subtypes and local AFS statistics
source("./scripts/window_analysis/subtype_abundances_byD.R")
source("./scripts/window_analysis/obs_vs_expected.R")





