#!/usr/bin/env python

# import packages
import numpy as np
import pandas as pd
import os
import os.path
import sys
from pathlib import Path
from scipy import stats

# retrieve command line arguments
args = sys.argv[1:]

# Access the arguments
sub = args[0] if len(args) > 0 else None
ses = args[1] if len(args) > 1 else None

# set directories
base_dir = "/data/CamRat/WHS_preprocessing/"
input_dir = base_dir + "outputs/MIND_files/MTR_scaled/cerebral_cortex/input/"
output_dir = base_dir + "outputs/MIND_files/MTR_scaled/cerebral_cortex/output/"

# create output directory if it doesn't exist
path = Path(output_dir)
path.mkdir(parents = True, exist_ok = True)

# path to MIND functions
mind_path = base_dir + 'MIND'
sys.path.insert(1, mind_path)
from MIND import compute_MIND # make sure nibabel and scipy are installed
from MIND_helpers import get_KL, calculate_mind_network

# Ignore warnings
#import warnings
#warning.filterwarnings("ignore")

# set file
file =  sub + "_" + ses + ".csv"

print("Calculating MIND network for " + sub + " " + ses)

# calculate MIND network if it has not already been generated
if os.path.isfile(output_dir + file):
    print("MIND network has already been generated")
else:
    print("generating MIND network")
    
    # read csv
    input_file = input_dir + file
    voxel_df = pd.read_csv(input_file)
    
    # get list of regions
    regions = voxel_df['Label'].unique()
    
    # NEW: create empty df with columns as ROIs and number of rows equal to number of resamples
    n_resamples = 5000
    resampled_dataset = pd.DataFrame(np.zeros((n_resamples, len(regions))), columns = regions)
    
    # NEW: For each ROI, estimate distribution and fill in number of resamples
    for name, data in voxel_df.groupby('Label'):
        resampled_dataset[name] = stats.gaussian_kde(data['MTR']).resample(n_resamples)[0]
    
    # NEW: reformat resampled dataset for MIND calculation
    resampled_dataset = resampled_dataset.melt(var_name = 'Label', value_name = 'MTR')
    
    # run MIND *on resampled dataset*
    rat_mind = calculate_mind_network(resampled_dataset, ['MTR'], regions)
    
    # save output
    output_path = output_dir + file
    rat_mind.to_csv(output_path, index = False)

print("finished!")
