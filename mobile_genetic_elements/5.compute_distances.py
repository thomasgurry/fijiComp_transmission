import sys
import os
import numpy as np
import pandas as pd
import pickle as pkl
import scipy
import scipy.stats
import statsmodels.sandbox.stats.multicomp
from numpy.random import RandomState
import multiprocessing as mp


# Get list of genomes and load depths
sample_dict = pkl.load(open('sample_dict.pkl', 'rb'))  # dict containing sample metadata
depth_dict = pkl.load(open('compressed_depth_dict_stool.pkl', 'rb'))
full_list_of_genomes = depth_dict[depth_dict.keys()[0]].keys()

# Load pairs
network_type = 'households'
related_file = 'same_household_pairs.pkl'   # pairs of sampleIDs that are somehow related
related_pairs = pkl.load(open(related_file, 'rb'))
unrelated_file = 'different_household_pairs.pkl'  # pairs of sampleIDs that are unrelated
unrelated_pairs = pkl.load(open(unrelated_file, 'rb'))

# Get a list of pairs with data
good_related_pairs = []
good_unrelated_pairs = []

for pair in related_pairs:

    # Make sure they have data
    try:
        person1_stool_GID = sample_dict[pair[0]]['stool_GID']
        person2_stool_GID = sample_dict[pair[1]]['stool_GID']
        p1 = depth_dict[person1_stool_GID][full_list_of_genomes[0]]
        p2 = depth_dict[person2_stool_GID][full_list_of_genomes[0]]        
    except:
        continue

    good_related_pairs.append([person1_stool_GID, person2_stool_GID])
    
for pair in unrelated_pairs:

    # Make sure they have data
    try:
        person1_stool_GID = sample_dict[pair[0]]['stool_GID']
        person2_stool_GID = sample_dict[pair[1]]['stool_GID']
        p1 = depth_dict[person1_stool_GID][full_list_of_genomes[0]]
        p2 = depth_dict[person2_stool_GID][full_list_of_genomes[0]]        
    except:
        continue

    good_unrelated_pairs.append([person1_stool_GID, person2_stool_GID])


print "Data loaded.  Launching MP module..."


# Start iterations
N_iterations = 100
related_counts = []
unrelated_counts = []
for N in range(N_iterations):

    # Downsample non-related pairs and find good set
    np.random.shuffle(unrelated_pairs)
    good_unrelated_pairs = []

    for pair in unrelated_pairs:
        try:
            person1_stool_GID = sample_dict[pair[0]]['stool_GID']
            person2_stool_GID = sample_dict[pair[1]]['stool_GID']
            p1 = depth_dict[person1_stool_GID][genome]
            p2 = depth_dict[person2_stool_GID][genome]
            good_unrelated_pairs.append([person1_stool_GID, person2_stool_GID])
        except:
            continue

        if len(good_unrelated_pairs) == len(good_related_pairs):
            break

    genome_counts = {}

    counter = 0
    for genome in full_list_of_genomes:

        genome_counts_related = 0
        genome_counts_unrelated = 0
        
        related_pairs_counted = 0
        # Get distances for related pairs
        for pair in good_related_pairs:

            p1 = depth_dict[pair[0]][genome]
            p2 = depth_dict[pair[1]][genome]

            Outlier1_indices = (p1 < np.median(p1)/1000) & (np.array(p2) >= np.median(p2))
            Outlier1 = np.shape(np.where(Outlier1_indices == True))[1]
            Outlier2_indices = (p2 < np.median(p2)/1000) & (np.array(p1) >= np.median(p1))
            Outlier2 = np.shape(np.where(Outlier2_indices == True))[1]
            diff = (Outlier1 + Outlier2)/float(len(p1))

            if diff == 0 and np.median(p1)>=10 and np.median(p2)>=10:
                genome_counts_related += 1
            related_pairs_counted += 1

            
        unrelated_pairs_counted = 0
        # Get distances for unrelated pairs
        for pair in good_unrelated_pairs:

            p1 = depth_dict[pair[0]][genome]
            p2 = depth_dict[pair[1]][genome]

            Outlier1_indices = (p1 < np.median(p1)/1000) & (np.array(p2) >= np.median(p2))
            Outlier1 = np.shape(np.where(Outlier1_indices == True))[1]
            Outlier2_indices = (p2 < np.median(p2)/1000) & (np.array(p1) >= np.median(p1))
            Outlier2 = np.shape(np.where(Outlier2_indices == True))[1]
            diff = (Outlier1 + Outlier2)/float(len(p1))

            if diff == 0 and np.median(p1)>=10 and np.median(p2)>=10:
                genome_counts_unrelated += 1
            unrelated_pairs_counted += 1

            if unrelated_pairs_counted == related_pairs_counted:
                break

        genome_counts[genome] = {}
        genome_counts[genome]['related'] = genome_counts_related
        genome_counts[genome]['unrelated'] = genome_counts_unrelated
        
        counter += 1

    # Count genomes in favour of related VS unrelated
    related_genomes = 0
    unrelated_genomes = 0

    for genome in genome_counts:
        if genome_counts[genome]['related'] > genome_counts[genome]['unrelated']:
            related_genomes += 1
        elif genome_counts[genome]['related'] < genome_counts[genome]['unrelated']:
            unrelated_genomes += 1
        else:
            continue

    related_counts.append(related_genomes)
    unrelated_counts.append(unrelated_genomes)


x = {'related_counts': related_counts, 'unrelated_counts': unrelated_counts}
pkl.dump(x, open('summary_stats_' + network_type + '.pkl', 'wb'))
