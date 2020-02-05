#! /usr/bin/env python

import msprime, pyslim
import numpy as np
import pandas as pd
import re
import sys


infile = sys.argv[1]
outfile = re.search("(.*).trees", infile).group(1)

ts = pyslim.load(infile).simplify()
samples = ts.num_samples

breaks = np.zeros(ts.num_trees + 1)
p1_hom = np.zeros(ts.num_trees + 1)
p2_hom = np.zeros(ts.num_trees + 1)

def getHomozygosity(tree, child, root_range):
	count = 0
	hom_count = 0
	for root in root_range:
		if tree.is_descendant(child, root):
			count += 1
		if tree.is_descendant(child+1, root):
			count += 1
		if count == 2:
			hom_count = 1
			break
	return hom_count

def homProp(child_range, root_range, tree):
	hom_count = 0
	for child in child_range:
		hom_count += getHomozygosity(tree, child, root_range)
	child_end = child_range[-2]+2
	child_start = child_range[0]
	hom_prop = hom_count / ((child_end - child_start) / 2)
	return hom_prop

p1_range = range(0,4)
p2_range = range(4,8)
child_range = range(8, samples, 2)

for tree in ts.trees(sample_counts=True):
	p1_prop = homProp(child_range, p1_range, tree)
	p2_prop = homProp(child_range, p2_range, tree)
	breaks[tree.index] = tree.interval[0]
	p1_hom[tree.index] = p1_prop
	p2_hom[tree.index] = p2_prop
	
breaks[-1] = ts.sequence_length
p1_hom[-1] = p1_hom[-2]
p2_hom[-1] = p2_hom[-2]

hom_df = pd.DataFrame(np.stack((breaks, p1_hom, p2_hom), axis = -1), columns = ["GenomicPosition", "P1HomProportion", "P2HomProportion"])

outname = f"{outfile}_homozygosityproportions.csv"

hom_df.to_csv(outname)