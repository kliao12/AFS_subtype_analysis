# Analysis for Effect of Mutation Subtypes on the Allele Frequency Spectrum and Population Genetics Inference # 

Step 0) Download zipped analysis directory from: zenodo.org. Use tar -xvzf afs_analysis.tar to unzip file. Major files in directories include: 
a) ./data
  - ./data/mst_sfs/: Text file for each of 96 subtypes 
  - ./gw_sites_info.txt: Text file containing each variant info used to make AFS
  - ./all_MST_100kb_counts.txt: Text file containing data for 100Kb window analysis (D, mst counts, etc) 
  - ./data/DaDi/bootstrap_SFS: 100 bootstrapped AFS for each subtype 
  - ./data/simulated_neutral/all_*.txt: Simulated neutral AFS using theta estimates for A_C.AAA and C_T.ACG
  
b) ./output
  - sigma_mat.txt: Large sigma matrix used for covariance calculations in D-2 statistic 
  - ./output/D2_stat/by_subtype: Output of D-2 statistic for each subtypes AFS
  - ./output/DaDi: Output from running dadi 10 times on each subtypes AFS.
  
c) ./scripts 
  - entire_analysis.R: R script containing entire pipeline 

Step 1) Software Needed to run all analysis completely from scratch.  
- R libraries needed can be found in ./scripts/libraries_needed.R
- DaDi (https://dadi.readthedocs.io/en/latest/user-guide/installation/)
- SLURM like parallel computing envrionment 




