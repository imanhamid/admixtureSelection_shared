#! /usr/bin/env python

'''
Calculates proportion of individuals (NOT chromosomes) with
p2 global ancestry > 0.75. Commented code will output txt file with
list of global ancestry proportions for all individuals.

This script is specific to the output from CV-admixture_chr1-rates.slim
and would need to be slightly modified to work for other chromosomes or population sizes.

Takes trees file as input, and prints to stdout the proportion. This is to make it easy
to run this script on many trees files in parallel, and then cat the output to one file
after all jobs are complete.

Usage:
globalancestry.py filename.trees
'''

import msprime, pyslim
import numpy as np
import sys
import re
##import json

infile = sys.argv[1]

#these are commented out because they are specific to how I was naming my files
#selection coefficient and seed number are not necessary to run this code
##s = re.search("seed-.*_s-(.*)\.trees", infile).group(1)
##seed = re.search("seed-(.*)_s-.*\.trees", infile).group(1)

ts = pyslim.load(infile).simplify()

#size of chr1
chrom_length = 249904549.

#table collection that makes up the full tree sequence
tablecoll=ts.dump_tables()

#map p1 ancestors to all samples (chromosomes) in admixed population
#outputs an EdgeTable with
#left genomic position, right genomic position, ancestor ID, and sample ID
ancestor_link = tablecoll.map_ancestors(range(20000,40000), range(10000))

#calculate interval over which each individual has ancestry from
#a specific ancestor
interval = ancestor_link.right - ancestor_link.left

##global_list = []
prop_count=0

#iterate over all individuals
#each individual has two chromosomes, so step size is 2
for child in range(20000,40000,2):
	#indices from EdgeTable for each individual (2 per individual)
	child_indx = np.argwhere(np.logical_or(ancestor_link.child==child, ancestor_link.child==child+1))
	
	#get global ancestry proportion:
	#sum over intervals corresponding to p1 ancestry and
	#divide by chrom_length*2 (because two chromosomes per individual)
	global_ancestry = interval[child_indx].sum()/(chrom_length * 2)
	
##	global_list.append(global_ancestry)
	
	#count proportion of individuals in the population
	#with p1 global ancestry > 0.75
	if global_ancestry > 0.75:
		prop_count+=1

prop = prop_count/10000

#print proportion of individuals with p1 global ancestry > 0.75 to stdout
##print(f"{s}\t{seed}\t{prop}")
print(f"{prop}")

##with open(f"./seed-{seed}_s-{s}_globalancestry.txt", "w") as f:
##	f.write(json.dumps(global_list))
