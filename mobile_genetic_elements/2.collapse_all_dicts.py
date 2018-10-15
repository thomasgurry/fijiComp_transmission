import os, sys
import pickle as pkl
import numpy as np
import multiprocessing as mp

# Script to combine all depth dicts into a single megadict with all 1kb regions in it

# Get full list of scaffolds and segments
scaffold_dict = {}
depth_dict_files = os.listdir('depth_dicts/')
for fn in depth_dict_files:
    x = pkl.load(open('depth_dicts/' + fn, 'rb'))
    for scaffold in x:
        if scaffold not in scaffold_dict:
            scaffold_dict[scaffold] = []
        for segment in x[scaffold]:
            if segment not in scaffold_dict[scaffold]:
                scaffold_dict[scaffold].append(segment)

# Get full list of subjects
subjects = []
for fn in depth_dict_files:
    GID = fn.split('_')[0]
    subjects.append(GID)
            
# Make megadict of dict[subject][scaffold][segment] = median depth, with 0 depth if not present
megadict = {GID: {} for GID in subjects}
for scaffold in scaffold_dict:
    for GID in megadict.keys():
        megadict[GID][scaffold] = {}
    for segment in scaffold_dict[scaffold]:
        for GID in megadict.keys():
            megadict[GID][scaffold][segment] = 0

for fn in depth_dict_files:
    x = pkl.load(open('depth_dicts/' + fn, 'rb'))
    GID = fn.split('_')[0]
    for scaffold in x:
        for segment in x[scaffold]:
            megadict[GID][scaffold][segment] = x[scaffold][segment]
            
# Pickle megadict
pkl.dump(megadict, open('depths_all_subjects.pkl', 'wb'))
