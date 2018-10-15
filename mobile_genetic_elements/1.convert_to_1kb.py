import os, sys
import pickle as pkl
import numpy as np
import multiprocessing as mp


def collapse_depths(fn):

    GID = fn.split('_')[0]

    # Read samtools depth output
    with open('depth_files_stool_f95_10kb/' + fn, 'r') as fid:
        all_lines = fid.readlines()
    
    # For each scaffold, create 1kb segments starting from 0 and count the median depth
    scaffold_dict = {}
    counter = 0
    for line in all_lines:
        spl = line.split()
        scaffold = spl[0]
        position = int(spl[1])
        segment_start = 1000 * (position/1000) # integer division
        depth = int(spl[2])
        if scaffold not in scaffold_dict:
            scaffold_dict[scaffold] = {}
        if segment_start not in scaffold_dict[scaffold]:
            scaffold_dict[scaffold][segment_start] = []
        scaffold_dict[scaffold][segment_start].append(depth)
        #counter += 1
        #if counter % 1000000 == 0:
        #    print str(counter) + ' of ' + str(len(all_lines))

    # Extract median depth for each segment
    depth_dict = {}
    for scaffold in scaffold_dict:
        if scaffold not in depth_dict:
            depth_dict[scaffold] = {}
        for segment in scaffold_dict[scaffold]:
            depth_dict[scaffold][segment] = np.median(scaffold_dict[scaffold][segment])

    # Dump to pickle
    pkl.dump(depth_dict, open('depth_dicts/' + GID + '_1kb_regions.pkl', 'wb'))
    

# Run in parallel on each depth file
filenames = os.listdir('depth_files_stool_f95_10kb/')
cpu_count = mp.cpu_count()
manager = mp.Manager()
#pool = mp.Pool(cpu_count)
pool = mp.Pool(2)

for fn in filenames:
    print fn
    
    pool.apply_async(collapse_depths, (fn,))

pool.close()
pool.join()
