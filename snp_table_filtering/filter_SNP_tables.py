import cPickle
import numpy as np
from scipy.stats import poisson
from scipy.stats import chisquare
import os
import pandas as pd

def load_sample_names():
    samples = np.array([line.rstrip() for line in open('sampleIDs.txt')])
    return samples


def load_kpileup(fn):
    # load the kpileup alignment
    # x[gene] = alignment (numpy array)
    x = cPickle.load(open(fn))
    return x


def remove_100bp(x):
    # remove 100 bp from beg/end of each gene
    # and organise them into separate arrays of size
    # MxNxK for M subjects, N positions, K nucleotides
    z = {}
    genes = sorted(x.keys())
    m = np.shape(x.values()[0])[0]
    k = 4
    for gene in genes:
        n = np.shape(x[gene])[1]-200
        if n < 0:
            continue
        y = np.zeros([m, n, k])
        y[:, :, :] = x[gene][:, 100:-100, :]
        z[gene] = y    
    return z

def concat_genes(x):
    # concat gene alignments for each subject
    # remove 100 bp from beg/end of each gene
    # returns MxNxK alignment (numpy array)
    genes = sorted(x.keys())
    m = np.shape(x.values()[0])[0]
    n = sum([(np.shape(x[gene])[1]) for gene in genes])
    k = 4
    y = np.zeros([m, n, k])
    beg = 0
    end = 0
    for gene in genes:
        end += (np.shape(x[gene])[1])
        y[:,beg:end,:] = x[gene][:,:,:]
        beg = end
    return y


def get_snps(x):
    # return an alignment containing only snps
    pos =  ((x > 0).sum(axis=2) > 1).sum(axis=0) > 0
    if sum(pos) > 0:
        return x[:, pos, :], pos
    else:
        return '', ''


def coverage(x):
    # calculate coverage for each subject
    # returns MxN numpy array
    return x.sum(axis=2)


def z_coverage(x):
    # calculate standardized coverage for each subject
    cov = coverage(x)
    zcov = (cov - cov.mean(axis=1)[:, np.newaxis]) / cov.std(axis=1)[:, np.newaxis]
    return zcov


def get_nonzero_subjects(x):
    # get all subjects that have full alignments
    # returns Mx1 numpy array
    cov = coverage(x)
    return ((cov > 0).sum(axis=1) > 0)


def remove_atypical_coverage(x):
    # calculate standardized coverage
    # remove sites that are >1 std from mean
    zcov = z_coverage(x)
    x[((zcov < -1.5) | (1.5 < zcov)),:] = 0
    return x


def remove_zero_subjects(x):
    mi = get_nonzero_subjects(x)
    x = x[mi,:,:]
    return x, mi


def remove_nan_columns(x):
    # remove columns that have only NAs
    cov = coverage(x)
    ni = (cov > 0).sum(axis=0) > 0
    x = x[:,ni,:]
    return x, ni


def floor_zeros(x, fzeros=0.30):
    # if a person has (>fzeros) positions with 0 coverage, remove their alignment from x
    # for each subject, calculate number of sites with zero coverage
    cov = coverage(x)
    subjects = ((cov == 0).sum(axis=1) > fzeros*np.shape(x)[1])
    x[subjects, :, :] = 0
    return x


def filter_multicopy_genes(x):
    # Removes genes that appear several times.
    # Takes as input the counts dict from 'load_kpileup'.
    # Returns a filtered array

    if x.keys()[0].split(';')[0] == 'ST' or x.keys()[0].split(';')[0] == 'SA':
        name_index = 1
    else:
        name_index = 0
    
    # Look for genes that appear more than once
    gene_dict = {}
    for gene in x:
        gene_name = gene.split(';')[name_index]
        if gene_name not in gene_dict:
            gene_dict[gene_name] = []
        gene_dict[gene_name].append(gene)
    singlecopy_genes = []
    multicopy_genes = []
    for gene_name in gene_dict:
        if len(gene_dict[gene_name]) > 1:
            multicopy_genes.append(gene_name)
        else:
            singlecopy_genes.append(gene_name)

    # If more than two genes are present in multiple copies, exit routine and
    # ignore this genome
    if len(multicopy_genes) > 2:
        return None

    # Write single copy genes back to new dict and return
    y = {}
    for gene in x:
        gene_name = gene.split(';')[name_index]
        if gene_name in singlecopy_genes:
            y[gene] = x[gene]

    return y
            

def check_gene_coverage_distribution(x):
    # Calculates chi square goodness of fit from sequencing depth data (counts per locus)
    # for each gene, compared to a Poisson distribution of the same mean.
    # The purpose of this test is to test whether read recruitment is Poisson or whether
    # there is underlying structure in a region of the gene that is recruiting reads from
    # other species, which could affect SNP calling.

    nsubjects = np.shape(x[x.keys()[0]])[0]

    # Loop through each gene
    chisq_pvalues = {gene: [] for gene in x}
    for gene in x:
        counts = {i: np.array([]) for i in range(nsubjects)}
        gene_counts = x[gene]

        # Compute background distribution for single copy genes for each subject
        for i in range(nsubjects):
            counts[i] = np.sum(gene_counts[i], axis=1)
            
            if np.nonzero(counts[i])[0].size > 0:
                hist, bin_edges = np.histogram(counts[i], bins=range(int(counts[i].max())+2)) # if counts observed = 0 or 1 only, want bin edges to be [0 1 2]
                poisson_mean = np.dot(hist, bin_edges[:-1])/float(np.sum(hist))

                # Compute empirical distribution
                empirical_freq = np.divide(hist, float(np.sum(hist)))
                empirical_freq = np.append(empirical_freq, np.array([0.0])) # add zero counts for the last bin so there are at least 3 bins for the chisquare test

                # Compute equivalent expected distribution
                poisson_pmf = poisson.pmf(range(bin_edges.max() + 20), poisson_mean)                    
                expected_freq = poisson_pmf[:len(bin_edges)]
                expected_freq[-1] += np.sum(poisson_pmf[len(bin_edges):])

                # Compute chi square goodness of fit - 1 extra degree of freedom because computed Poisson mean from data
                chisq, pval = chisquare(empirical_freq, f_exp=expected_freq, ddof=1)
                chisq_pvalues[gene].append(pval)
                
    # Return median p-values
    median_pvalues = {gene: np.median(chisq_pvalues[gene]) for gene in x}
    return median_pvalues    


def filter_genes_by_coverage_distribution(x):
    # Calculates chi square goodness of fit from sequencing depth data (counts per locus)
    # for each gene, compared to a Poisson distribution of the same mean.
    # Returns a filtered dict with genes with mean p-value cut-off lower than 0.05
    # thrown out.
    # The purpose of this test is to test whether read recruitment is Poisson or whether
    # there is underlying structure in a region of the gene that is recruiting reads from
    # other species, which could affect SNP calling.

    nsubjects = np.shape(x[x.keys()[0]])[0]

    # Loop through each gene
    chisq_pvalues = {gene: [] for gene in x}
    for gene in x:
        counts = {i: np.array([]) for i in range(nsubjects)}
        gene_counts = x[gene]

        # Compute background distribution for single copy genes for each subject
        for i in range(nsubjects):
            counts[i] = np.sum(gene_counts[i], axis=1)
            
            if np.nonzero(counts[i])[0].size > 0:
                hist, bin_edges = np.histogram(counts[i], bins=range(int(counts[i].max())+2)) # if counts observed = 0 or 1 only, want bin edges to be [0 1 2]
                poisson_mean = np.dot(hist, bin_edges[:-1])/float(np.sum(hist))

                # Compute empirical distribution
                empirical_freq = np.divide(hist, float(np.sum(hist)))
                empirical_freq = np.append(empirical_freq, np.array([0.0])) # add zero counts for the last bin so there are at least 3 bins for the chisquare test

                # Compute equivalent expected distribution
                poisson_pmf = poisson.pmf(range(bin_edges.max() + 20), poisson_mean)                    
                expected_freq = poisson_pmf[:len(bin_edges)]
                expected_freq[-1] += np.sum(poisson_pmf[len(bin_edges):])

                # Compute chi square goodness of fit - 1 extra degree of freedom because computed Poisson mean from data
                chisq, pval = chisquare(empirical_freq, f_exp=expected_freq, ddof=1)
                chisq_pvalues[gene].append(pval)
                
    # Compute median p-values
    median_pvalues = {gene: np.median(chisq_pvalues[gene]) for gene in x}

    # Remove genes with median p-values less than 0.05
    y = {}
    for gene in median_pvalues:
        if median_pvalues[gene] >= 0.05:
            y[gene] = x[gene]

    return y

        
def process_kpileup(fn, subjectcounts_df=None, bacterial_genes=None, out=''):
    # process a kpileup alignment
    x = load_kpileup(fn)
    number_of_genes_at_start = len(x.keys())
    M = load_sample_names()
    x = remove_100bp(x)
    if x == None or x == {}:
        print "Nothing left after removing 100bp on either end of each gene."
        return np.array([]), [], 0, 0, subjectcounts_df, []
    x = filter_multicopy_genes(x)
    if x == None  or x == {}:
        print "Nothing left after filtering multicopy genes."
        return np.array([]), [], 0, 0, subjectcounts_df, []
    x = filter_genes_by_coverage_distribution(x)
    if x == None or x == {}:
        print "Nothing left after filtering genes with non-Poisson coverage."
        return np.array([]), [], 0, 0, subjectcounts_df, []
    if x.keys()[0].split(';')[0] == 'ST' or x.keys()[0].split(';')[0] == 'SA':
        name_index = 1
    else:
        name_index = 0
    # Remove non-bacterial genes if provided, else assume archaea and accept all genes
    if bacterial_genes != None:
        genes_to_delete = []
        for gene in x:
            gene_name = gene.split(';')[name_index]
            if gene_name not in bacterial_genes:
                genes_to_delete.append(gene)
        for gene in genes_to_delete:
            del x[gene]
            
    # Count genes and store them
    try:
        genes_after_filtering = x.keys()
        number_of_genes_after_filtering = len(x.keys())
    except:
        pass
        
    # Return nothing if there are more than two genes present in multiple copies
    if x == None or x == {}:
        print "More than two genes with multiple copies detected.  Aborting."
        return np.array([]), [], 0, 0, subjectcounts_df, []

    # Add entry to subject counts DataFrame (number of subjects per gene, after filtering)
    if isinstance(subjectcounts_df, pd.DataFrame):
        column_names = subjectcounts_df.columns.tolist()
        newentry = pd.DataFrame([len(column_names)*[np.nan]], columns = column_names)
        genome = fn.split('/')[-1].split('.')[0]
        identity = fn.split('/')[-1].split('.')[1].split('pid')[0]
        newentry['Genome'] = genome
        newentry['Alignment identity'] = identity
        for gene in x:
            gene_name = gene.split(';')[name_index]
            nsubjects = np.shape(x[gene])[0]
            newentry[gene_name] = 0
            for k in range(nsubjects):
                if np.sum(x[gene][k]) > 0:
                    newentry[gene_name] += 1
            
        subjectcounts_df = pd.concat([newentry, subjectcounts_df], axis=0)
        
    x = concat_genes(x)
    x, ni = get_snps(x)
    if x == '':
        return np.array([]), [], 0, 0, subjectcounts_df, []
    x, mi = remove_zero_subjects(x)
    M = M[mi]    
    if x == '':
        return np.array([]), [], 0, 0, subjectcounts_df, []
    if np.shape(x)[0] > 0:
        x, ni = remove_nan_columns(x)    
    return x,M, number_of_genes_at_start, number_of_genes_after_filtering, subjectcounts_df, genes_after_filtering
    

def check_coverage(fn, outfile=''):
    # process a kpileup alignment
    x = load_kpileup(fn)
    M = load_sample_names()
    x = remove_100bp(x)
    x = filter_multicopy_genes(x)
    median_pvalues = check_gene_coverage_distribution(x)
    return median_pvalues


import sys

identity=sys.argv[1]
if identity not in [90, 95, 97, 99]:
    exit

#identity=90

alignment_files = os.listdir('alignments/')
output_directory = 'output/'

# Summary dataframes initialization
columns_genecounts = ['Genome', 'Identity', 'Number of genes at start', 'Number of genes after filtering']
genecounts_df = pd.DataFrame(columns=columns_genecounts)
bacterial_genes = ['rplA', 'rplB', 'rplC', 'rplD', 'rplE', 'rplF', 'rplK', 'rplL', 'rplM',
                   'rplN', 'rplP', 'rplS', 'rplT', 'rpmA', 'rpsB', 'rpsC', 'rpsE', 'rpsI',
                   'rpsJ', 'rpsK', 'rpsM', 'rpsS', 'tsf', 'frr', 'infC', 'smpB', 'rpoB',
                   'nusA', 'dnaG', 'pgk', 'pyrG']
columns_subjectcounts = ['Genome', 'Alignment identity'] + bacterial_genes
subjectcounts_df = pd.DataFrame(columns=columns_subjectcounts)
genes_per_genome = {} # Track genes left after filtering for each genome

file_counter = 0
for fn in alignment_files:
    print fn

    # Check if the genome is archaeal
    if fn.split('.')[0] == '3054Methanobrevibacter' or fn.split('.')[0] == '3199Methanobrevibacter':
        x, M, number_of_genes_at_start, number_of_genes_after_filtering, subjectcounts_df, genes_after_filtering = process_kpileup('alignments/'+ fn, subjectcounts_df)
    else:
        x, M, number_of_genes_at_start, number_of_genes_after_filtering, subjectcounts_df, genes_after_filtering = process_kpileup('alignments/' + fn, subjectcounts_df, bacterial_genes)

    genome_name = fn.split('.pickle')[0]
    genome = fn.split('.')[0]
    genes_per_genome[genome_name] = genes_after_filtering
    newentry = pd.DataFrame([len(columns_genecounts)*[np.nan]], columns=columns_genecounts)
    newentry['Genome'] = genome
    newentry['Identity'] = fn.split('.')[1].split('pid')[0]
    newentry['Number of genes at start'] = number_of_genes_at_start
    newentry['Number of genes after filtering'] = number_of_genes_after_filtering
    genecounts_df = pd.concat([newentry, genecounts_df], axis=0)
    
    if np.shape(x)[0] > 0:
        file_name = fn.split('.pickle')[0]
        cPickle.dump(x, open(output_directory + file_name + '.aln.pickle', 'w'))
        cPickle.dump(M, open(output_directory + file_name + '.ids.pickle', 'w'))

    file_counter += 1
    print str(file_counter) + " out of " + str(len(alignment_files)) + " processed."
    
