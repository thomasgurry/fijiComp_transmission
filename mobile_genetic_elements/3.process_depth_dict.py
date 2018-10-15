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
pkl.dump(genomeIDs, open('genome_scaffold_dict_stool.pkl', 'wb'))
print "Extracted all genome and scaffold IDs from dataset."

# Combine genome scaffolds into single arrays (72: [scaffold001, scaffold002, etc.])
# e.g. 72: [00002232322310000000]
#           scaffold1  scaffold2

final_dict = {GID: {} for GID in megadict.keys()}
segment_key = {genomeID: [] for genomeID in genomeIDs}
for GID in final_dict:
    print GID
    for genomeID in genomeIDs:
        depth_array = []
        segmentID_array = []
        for scaffoldID in genomeIDs[genomeID]:
            for segment in megadict[GID][genomeID + '_scaffold' + scaffoldID]:
                depth_array.append(megadict[GID][genomeID + '_scaffold' + scaffoldID][segment])
                segmentID_array.append(genomeID + '_scaffold' + scaffoldID + '_segment' + str(segment))
        final_dict[GID][genomeID] = depth_array
        segment_key[genomeID] = segmentID_array

pkl.dump(segment_key, open('segment_map_for_compressed_depth_dict.pkl', 'wb'))
pkl.dump(final_dict, open('compressed_depth_dict_stool.pkl', 'wb'))
