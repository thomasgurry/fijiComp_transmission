"""

Script to calculate Manhattan distances at the level of SNPs 
between pairs of subjects for a given set of bacteria.

"""

import os
import numpy as np
import pickle as pkl
import multiprocessing as mp


def make_presence_absence(snp_table):
    # Takes as input a SNP array [Nx4] of allele counts [A, G, C, T]
    # and returns the same array with just 1's or 0's (presence/absence).
    # Requires at least 2 counts of an allele to return a 'presence' count.
    presence_absence_array = np.zeros(np.shape(snp_table))
    inds = snp_table > 0
    presence_absence_array[inds] = 1
    return presence_absence_array


def make_presence_absence_dominant_allele_only(snp_table):
    # Takes as input a SNP array [Nx4] of allele counts [A, G, C, T]
    # and returns the same array with just 1's or 0's (presence/absence).
    # Only the dominant allele is returned.
    presence_absence_table = np.zeros(np.shape(snp_table))
    for i in range(np.shape(snp_table)[0]):
        if np.sum(snp_table[i, :]) > 0:
            presence_absence_table[i, np.argmax(snp_table[i, :])] = 1            
    return presence_absence_table
    
def manhattan_distance(snp_table_1, snp_table_2):
    # Calculates Manhattan distance of presence/absence at each SNP,
    # then normalizes by number of SNPs to return an average distance
    # per locus for that genome.

    # Presence absence allele vector for each subject
    subject1 = make_presence_absence_dominant_allele_only(snp_table_1)
    subject2 = make_presence_absence_dominant_allele_only(snp_table_2)
    loci_to_keep = np.all(np.array([np.sum(subject1, axis=1) > 0, np.sum(subject2, axis=1) > 0]), axis=0)
    dist = np.sum(np.abs(subject1[loci_to_keep, :] - subject2[loci_to_keep, :]))
    return dist/len(loci_to_keep)


def parse_genome(manhattan_distances_dict, genome, identity):
    # Takes as input a bacterial genome name (e.g. 'BACT_1'), identity level (e.g. 95) and list of sampleIDs (in the order
    # in which they are present in the resulting distance matrix).
    # Computes Manhattan distances between each subject in the associated alignment file.

    # Load SNP tables and associated sample IDs
    alignment_file = '/path/to/filtered/SNP_tables/' + genome + '.' + str(identity) + 'pid.aln.pickle'
    sampleID_file = '/path/to/filtered/SNP_tables/' + genome + '.' + str(identity) + 'pid.ids.pickle'
    snp_table = pkl.load(open(alignment_file, 'rb'))
    sampleID_list = pkl.load(open(sampleID_file, 'rb'))

    # Read in full list of N sample IDs and create a matrix of NxN distances for each genome
    with open('sampleIDs.txt', 'r') as fid:
        al = fid.readlines()
    sampleIDs = [line.rstrip('\n') for line in al]    
    tmp_array = np.zeros((len(sampleIDs), len(sampleIDs)))
    tmp_array[:,:] = np.NAN
    #manhattan_distances_dict[genome] = np.zeros((len(sampleIDs), len(sampleIDs)))
    #manhattan_distances_dict[genome][:,:] = np.NAN

    # Figure out indices of samples in full sample list to know correct indices to use in distance matrices
    sample_indices = [sampleIDs.index(sampleID) for sampleID in sampleID_list]
    
    # Compute pairwise Manhattan distances for each subject pair
    for i in range(len(sample_indices)):

        snps_i = snp_table[i]

        for j in range(len(sample_indices)):

            snps_j = snp_table[j]

            # Compute Manhattan distance
            dist = manhattan_distance(snps_i, snps_j)
            tmp_array[sample_indices[i], sample_indices[j]] = dist

    manhattan_distances_dict[genome] = tmp_array
    
# Step 0 - obtain full set of bacterial genome alignments
all_alignments = os.listdir('/path/to/filtered/SNP_tables/')
alignments_90 = [alignment_file for alignment_file in all_alignments if alignment_file.split('.')[1] == '90pid' and alignment_file.split('.')[2] == 'aln']
genomes = [alignment_file.split('.')[0] for alignment_file in alignments_90]

# Step 1 - read in full list of N sample IDs and create a matrices of NxN distances for each genome
with open('sampleIDs.txt', 'r') as fid:
    al = fid.readlines()
sampleIDs = [line.rstrip('\n') for line in al]

# Step 2 - for each genome, compute pairwise Manhattan distance between subject
cpu_count = mp.cpu_count()
manager = mp.Manager()

identity='90'
manhattan_distances = manager.dict()
pool = mp.Pool(cpu_count)

for genome in genomes[identity]:

    pool.apply_async(parse_genome, (manhattan_distances, genome, identity,))

pool.close()
pool.join()

# Write distances to file
pkl.dump(dict(manhattan_distances), open('manhattan_distances.pkl', 'wb'))
