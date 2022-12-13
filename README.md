# Analysis for Effect of Mutation Subtypes on the Allele Frequency Spectrum and Population Genetics Inference # 

Step 1) Download zipped analysis directory from: https://zenodo.org/record/7430597#.Y5fphi1h1pQ. Use tar -xvzf afs_analysis.tar to unzip file. Major files in directories include: 

a) ./data
  - ./data/mst_sfs/: Text file for each of 96 subtypes. File all_subtypes_afs.txt contains all 96 AFS in one file 
  - ./gw_sites_info.txt: Text file containing each variant info genome wide used in analysis 
  - ./all_MST_100kb_counts.txt: Text file containing data for 100Kb window analysis (D, mst counts, etc) 
  - ./data/DaDi/bootstrap_SFS: 100 bootstrapped AFS for each subtype 
  - ./data/simulated_neutral/all_*.txt: Simulated neutral AFS using theta estimates for A_C.AAA and C_T.ACG
  
b) ./output
  - sigma_mat.txt: Large sigma matrix used for covariance calculations in D-2 statistic 
  - ./output/D2_stat/by_subtype: Output of D-2 statistic for each subtypes AFS
  - ./output/DaDi: Output from running dadi 10 times on each subtypes AFS.
  
c) ./scripts 
  - entire_analysis.R: R script containing entire pipeline sourcing each relevant analysis file 
  - afs_paper_functions.R: R script with functions used

Step 2) Install software needed to run all analysis completely from scratch.  
- R libraries needed can be found in ./scripts/libraries_needed.R
- DaDi (https://dadi.readthedocs.io/en/latest/user-guide/installation/)
- SLURM like parallel computing envrionment 

Step 3) Run analysis by going through entire_analysis.R in sections.
- Create AFS for each of 96 mutation subtypes. Outputs separate txt file for each and a single with all 96 subtypes.
- Create dataframe of genome wide AFS statistics for each of 96 subtypes. Includes D, D-2, % singles, doubles, triples, etc 
- AFS heterogeneity across subtypes and signals of recurrent mutations and gene conversion shaping the AFS 
  - Signals of mutation rate/recurrent mutations on AFS. Look at singleton/doubleton ratio by mutation rate
  - Signals of gBGC on AFS. Split into WS, SW, and indifferent mutation types. Compare D and D-2 across subtypes  
- Effect of genome wide AFS heterogeneity across subtypes on demographic inference using DaDi
- Effect of subtype composition and local genomic factors on the local AFS 
 



