# Analysis for Effect of Mutation Subtypes on the Allele Frequency Spectrum and Population Genetics Inference 

Step 0) Download zipped analysis directory from: zenodo.org. Use tar -xvzf afs_analysis.tar to unzip file. Major files in directories include: 
- ./data
  - 
- ./output
  - sigma_mat.txt: Large sigma matrix used for covariance calculations in D-2 statistic 
  - ./output/D2_stat: 
  - ./output/DaDi: Output from running dadi 10 times on each subtypes AFS.
- ./scripts 
  - entire_analysis.R: R script containing entire pipeline 

Step 1) Software Needed to run all analysis completely from scratch.  
- R libraries needed can be found in ./scripts/libraries_needed.R
- DaDi (https://dadi.readthedocs.io/en/latest/user-guide/installation/)
- SLURM like parallel computing envrionment 




