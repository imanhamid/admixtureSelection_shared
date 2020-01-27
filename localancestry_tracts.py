#! /usr/bin/env python

'''
Calculates from tree sequence file tract lengths for chromosomes 
that have the beneficial allele at the DARC locus.
This script is specific to the output from CV-admixture_chr1-rates.slim
and would need to be slightly modified to work for other chromosomes or population sizes.

Takes trees file as input, outputs txt file listing tract lengths
and their start/end posiitons

Usage:
localancestry_tracts.py filename.trees
'''

import msprime, pyslim
import numpy as np
import re
import sys

infile = sys.argv[1]
chrom_length = 249904549.
outfile = re.search("(.*)\.trees", infile).group(1)

##position of the beneficial locus
bene_locus = 159174683.

##read in trees file, find tree at the beneficial locus
ts = pyslim.load(infile).simplify()
bene_tree = ts.at(bene_locus)

##get sample IDs that have p1 ancestry at beneficial locus
bene_samples = []

for root in bene_tree.roots:
	if bene_tree.population(root)==0:
		for sample in bene_tree.samples(root):
			if sample >= 20000:
				bene_samples.append(sample)

##table collection that makes up the full tree sequence
tablecoll = ts.dump_tables()
##mapping ancestors to samples that had the beneficial allele
bene_ancestors = tablecoll.map_ancestors(bene_samples, range(10000))

##get local ancestry tract lengths surrounding the beneficial locus
##for each sample that had the beneficial trait

#tract_arrays = []
#tract_lengths = []
with open(f"{outfile}_tractlengths.txt", 'w') as out:
	out.write("tractLength\tstart\tend\n")
	for child in bene_samples:
		child_indx = np.argwhere(bene_ancestors.child==child) ##indices from EdgeTable for each specific sample
		left = []
		right = []
		
		##pull out start and end coords for each edge for this sample	
		for indx in child_indx:
			left.append(bene_ancestors[int(indx)].left)
			right.append(bene_ancestors[int(indx)].right)
		
		left = np.array(left)
		right = np.array(right)
		
		##sort in ascending order	
		left.sort(axis=0)
		right.sort(axis=0)
		
		##get index for interval that includes the beneficial locus
		bene_indx = int(np.argwhere(np.logical_and(left<=bene_locus, bene_locus<=right)))
		
		##combine start and end coords
		child_trees = np.stack((left, right), axis=-1)
		del left, right ##i don't know if this actually frees up any memory...

		start = child_trees[bene_indx][0]
		end = child_trees[bene_indx][1]
		
		##move up forward one interval and check if interval is contiguous with ancestry tract from beneficial locus
		##if so, save new end coordinate for local ancestry tract length
		for i in range(bene_indx, len(child_trees)-1):
			if child_trees[i][1] == child_trees[i+1][0]:
				end = child_trees[i+1][1]
			else:
				break
		
		##same as before, but moving in reverse
		for i in reversed(range(0,bene_indx+1)):
			if child_trees[i][0] == child_trees[i-1][1]:
				start = child_trees[i-1][0]
			else:
				break
		
		tract_length = end - start
		
		##save tract lengths to list
	#	tract_arrays.append([child, start, end])
	#	tract_lengths.append(tract_length)
		out.write(f"{tract_length}\t{start}\t{end}\n")


