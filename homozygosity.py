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

def getDescendents(tree, root_range):
	offspring = []
	for root in root_range:
		for child in tree.samples(root):
			offspring.append(child)
	return offspring

def homProp(offspring, child_start, child_end):
	hom_count = 0
	for child in range(child_start, child_end, 2):
		if child in offspring and child+1 in offspring:
			hom_count += 1
	hom_prop = hom_count / ((child_end - child_start) / 2)
	return hom_prop

for tree in ts.trees(sample_counts=True):
	p1_range = range(0,4)
	p2_range = range(4,8)
	p1_prop = homProp(getDescendents(tree, p1_range), 8, samples)
	p2_prop = homProp(getDescendents(tree, p2_range), 8, samples)
	breaks[tree.index] = tree.interval[0]
	p1_hom[tree.index] = p1_prop
	p2_hom[tree.index] = p2_prop
breaks[-1] = ts.sequence_length
p1_hom[-1] = p1_hom[-2]
p2_hom[-1] = p2_hom[-2]

hom_df = pd.DataFrame(np.stack((breaks, p1_hom, p2_hom), axis = -1), columns = ["GenomicPosition", "P1HomProportion", "P2HomProportion"])

outname = f"{outfile}_homozygosityproportions.csv"

hom_df.to_csv(outname)