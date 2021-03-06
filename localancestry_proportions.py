#! /usr/bin/env python

'''
Calculates population-wide local ancestry proportions across the genome.
This script is not specific to any chromosome or population size,
so it should be usable for any simulation that outputs
a trees file.

This script is VERY slow. It takes ~2 hrs per trees file when simulating
entire Chr1 with population size 10000. I tried using multiprocessing,
and that just ended up being slower. Not sure how to improve this.

Takes trees file as input, and outputs csv file listing genomic positions
and p2 (Euro for CV) ancestry proportion. If you want p1 (WAfrican for CV)
ancestry prop, you have to calculate that from p2 proportions.

Usage:
localancestry_proportions.py filename.trees
'''

import msprime, pyslim
import numpy as np
import pandas as pd
import re
import sys


infile = sys.argv[1]
outfile = re.search("(.*).trees", infile).group(1)

ts = pyslim.load(infile).simplify()

#create zero arrays for genomic positions
#and corresponding ancestry proportions
breaks = np.zeros(ts.num_trees + 1)
ancestry = np.zeros(ts.num_trees + 1)

#iterate over trees (1 unique geneology/genomic interval)
for tree in ts.trees(sample_counts=True):
	subpop_sum, subpop_weights = 0, 0

	#iterate over ancestors [roots]
	for root in tree.roots:
		#for each ancestor, count descendents
		leaves_count = tree.num_samples(root) - 1

		#get the pop of ancestor (p1=0, p2=1)
		#subpop_sum = num chromosomes for each genomic position w/ p1 ancestry
		subpop_sum += tree.population(root) * leaves_count

		#subpop_weights = num of chromosomes in sample
		subpop_weights += leaves_count

	#genomic position corresponds to left end of tree interval
	breaks[tree.index] = tree.interval[0]

	#ancestry proportion = num chroms with p1 ancestry / num chroms
	ancestry[tree.index] = subpop_sum / subpop_weights

#fill chromosome end positions
breaks[-1] = ts.sequence_length
ancestry[-1] = ancestry[-2]

#transform arrays to pandas dataframe
ancestry_arr = np.stack((breaks, ancestry), axis = -1)
ancestry_df = pd.DataFrame(ancestry_arr, columns = ["GenomicPosition", "P2AncestryProportion"])

outname = f"{outfile}_ancestryproportions.csv"

#save to csv file
ancestry_df.to_csv(outname)