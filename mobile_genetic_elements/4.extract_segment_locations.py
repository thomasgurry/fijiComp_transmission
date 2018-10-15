import os, sys
import pickle as pkl
import numpy as np
import multiprocessing as mp

# Combine depth dict into single arrays of median depths

# Load dict
megadict = pkl.load(open('depths_all_subjects_stool.pkl', 'rb'))

# Extract all genome IDs and scaffold IDs
genomeIDs = {}
for GID in megadict:
    print GID
    for scaffold in megadict[GID]:
        genomeID = scaffold.split('_')[0]
        scaffoldID = scaffold.split('_scaffold')[1]
        if genomeID not in genomeIDs:
            genomeIDs[genomeID] = []
        if scaffoldID not in genomeIDs[genomeID]:
            genomeIDs[genomeID].append(scaffoldID)

segment_key = {genomeID: [] for genomeID in genomeIDs}
GID = megadict.keys()[0]
for genomeID in genomeIDs:
    segmentID_array = []
    for scaffoldID in genomeIDs[genomeID]:
        for segment in megadict[GID][genomeID + '_scaffold' + scaffoldID]:
            segmentID_array.append(genomeID + '_scaffold' + scaffoldID + '_segment' + str(segment))
    segment_key[genomeID] = segmentID_array

pkl.dump(segment_key, open('segment_map_for_compressed_depth_dict.pkl', 'wb'))
