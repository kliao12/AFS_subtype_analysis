##### Title: Script to run DaDi and generate uncertanties
##### Date: 10/11/22

import dadi
import pandas as pd
import numpy
from multiprocessing import Process
import os
import subprocess
import sys
import nlopt
import glob 

##### Make sure to run this from DaDi_3 conda environment
run_num = sys.argv[1]
print(run_num)

MST = sys.argv[2]
print(MST)

#Import data
all_chr = pd.read_csv("/net/wonderland/home/ksliao/zoellner_research/scripts/output/genomewide_fullanno.vcf", delimiter="\t", header=None)
all_chr.columns = ["CHR", "REF","ALT", "AA", "AC_correct", "ANNO", "MT", "KMER", "MST"]

all_chr_MST = all_chr[all_chr['MST'] == MST]

#Create site frequency spectrum using frequency counts
item_list = range(7113);

all_chr_MST_afs = pd.DataFrame(all_chr_MST['AC_correct'].value_counts().reindex(item_list, fill_value = 0))
all_chr_MST_afs.columns = ['ac']
all_chr_MST_afs.sort_index(inplace=True)

#Create an array from dataframe ac column values
#all_chr_MST_afs = all_chr_MST_afs.head(10)
all_chr_MST_ac_array = all_chr_MST_afs['ac'].values

data_fs = dadi.Spectrum(all_chr_MST_ac_array)
ns = data_fs.sample_sizes

#Define demographic model
demo_model = dadi.Demographics1D.three_epoch

# If the data is unfolded (the ancestral allele was known), as the example data is,
# wrap the demographic model in a function that adds a parameter to estimate the rate of misidentification.
#demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)

# Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

#Initialize grid points and initial values. Parameter boudns too 
pts_l = [max(ns)+120, max(ns)+130, max(ns)+140]

params = [0.5, 20, 1, .2]

# Define boundaries of optimization.
upper_bounds = [1, 150, 10, 10]
lower_bounds = [1e-3, 1e-3, 1e-3, 1e-3]

#Run DaDi
p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds,
                              lower_bound=lower_bounds)

print("begin optimizing model")
popt, ll_model = dadi.Inference.opt(p0, data_fs, demo_model_ex, pts_l,
                             lower_bound=lower_bounds,
                             upper_bound=upper_bounds,
                             algorithm=nlopt.LN_BOBYQA,
                             maxeval=600, verbose=0)
print("finished optimizing model")
print(popt)

#Calculate theta 
model_fs = demo_model_ex(popt, ns, pts_l)
theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, data_fs)

# Import bootstrap SFS to get Godambe uncertainties
mst_boots = glob.glob("/net/wonderland/home/ksliao/zoellner_research/AFS_project/for_submission/data/DaDi/bootstrap_SFS/" + MST + "/*.txt" )

print("begin estimating uncertainties")
boot_sfs = [pd.read_csv(fid, header=None) for fid in mst_boots]

dadi_sfs = [dadi.Spectrum(sfs.iloc[:,0].to_numpy()) for sfs in boot_sfs]

uncerts_adj = dadi.Godambe.GIM_uncert(demo_model_ex, pts_l, dadi_sfs, popt, data_fs, eps=0.001, multinom=False)
    
print("finished estimating uncertainties")
print(uncerts_adj)

#Output to dataframe
param_df = pd.DataFrame([[ MST, ll_model, theta0, popt[0], popt[1], popt[2], popt[3], uncerts_adj[0], uncerts_adj[1], uncerts_adj[2], uncerts_adj[3] ]],columns=['MST', 'log_lik', 'theta', 'nuB','nuF','TB', 'TF', 'nuB_sd', 'nuF_sd', 'TB_sd', 'TF_sd'])
param_df.to_csv('/net/wonderland/home/ksliao/zoellner_research/AFS_project/for_submission/output/DaDi/' + MST + '_three_epoch' + run_num + '.txt', sep = '\t')

