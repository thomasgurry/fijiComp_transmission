# fijiComp_transmission
Scripts from Brito, Gurry et. al for computing transmission-related metrics from shotgun metagenomics data.

## SNPs from core (AMPHORA) genes
### SNP table filtering
Once alignments to core genes of reference genomes have been generated, they must be transformed into alignment tables for each reference genome, in the form of a Python dictionary containing individual core genes as its highest level keys, where for each core gene there is a M x N x 4 numpy array, for M subjects, N loci, and four different alleles (A, G, C and T).  Filtering steps described in the article to obtain SNP tables for downstream analyses are generated using the 'filter_SNP_tables.py' script.

### Distance calculations
Pairwise average Manhattan distances per SNP between individual subjects' SNP tables can be computed using the 'calculate_manhattan_distance_dominant_allele_parallel.py' script.

## Mobile genetic elements - computing summary statistics from presence/absence of 1kb mobile elements
Starting with outputs of 'samtools depth' (to get depth at each position), after modifying the appropriate paths in the scripts and organising metadata, run the following scripts in sequence:

**1.convert_to_1kb.py
**2.collapse_all_dicts.py
**3.process_depth_dict.py
**4.extract_segment_locations.py
**5.compute_distances.py
