#!/usr/bin/env python

# import packages
import numpy as np
import pandas as pd
import os
import os.path
import sys
from pathlib import Path

# retrieve command line arguments
args = sys.argv[1:]

# Access the arguments
sub = args[0] if len(args) > 0 else None
ses = args[1] if len(args) > 1 else None

# set directories
base_dir = "/data/CamRat/WHS_preprocessing/"
input_dir = base_dir + "outputs/MIND_files/MTR/GM_only_noBrainstemCerebellum/input/"
output_dir = base_dir + "outputs/MIND_files/MTR/GM_only_noBrainstemCerebellum/output/"

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

	# remove hemisphere from Label (turn off to split by hemisphere)
	voxel_df['Label'] = voxel_df['Label'].str.rsplit('_', n=1).str[0]

	# get list of regions
	regions = voxel_df['Label'].unique()

	# run MIND
	rat_mind = calculate_mind_network(voxel_df, ['MTR'], regions)

	# save output
	output_path = output_dir + file
	rat_mind.to_csv(output_path, index=False)

print("finished!")

