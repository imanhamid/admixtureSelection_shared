#! /usr/bin/env python

import msprime, pyslim
import numpy as np
import pandas as pd
import sys
import re

infile = sys.argv[1]
outfile = re.search("(.*).trees", infile).group(1)
ts = pyslim.load(infile).simplify()
chrom_length = ts.sequence_length
bene_locus = [site.position for site in ts.sites()][0]
samples = ts.num_samples
tablecoll = ts.dump_tables()
ancestor_link = tablecoll.map_ancestors(range(8,samples), range(4))
del tablecoll

def getChildArray(child, EdgeTable):
	indices = np.argwhere(EdgeTable.child==child)
	child_arr =[]
	for index in indices:
		child_arr.append([EdgeTable[int(index)].left, EdgeTable[int(index)].right])
		child_arr.sort(key = lambda x: x[0])
	return child_arr

def getTracts(child_arr):
	tracts = []
	start = -1
	max = -1
	for i in range(len(child_arr)):
		t = child_arr[i]
		if t[0] > max:
			if i != 0:
				tracts.append([start, max])
			max = t[1]
			start = t[0]
		else:
			if t[1] >= max:
				max = t[1]
	if max != -1 and [start,max] not in tracts:
		tracts.append([start,max])
	return tracts

def numSwitches(tracts, left, right):
	num_switches = 0
	if len(tracts) == 1:
		if tracts[0][0] == 0 and tracts[0][1] == chrom_length:
			num_switches = 0
		elif tracts[0][0] == 0 and tracts[0][1] != chrom_length:
			if left <= tracts[0][1] and tracts[0][1] <= right:
				num_switches+=1
		elif tracts[0][1] == chrom_length and tracts[0][0] != 0:
			if left <= tracts[0][0] and tracts [0][0] <= right:
				num_switches+=1
		elif left <= tracts[0][0] and tracts[0][1] <= right:
			num_switches = 2
	else:
		for i in range(len(tracts)):
			t = tracts[i]
			if left <= t[0] and t[1] <= right:
				if t[1] == chrom_length or t[0] == 0:
					num_switches += 1
				else:
					num_switches += 2
			elif left >= t[0] and t[1] >= right:
				num_switches += 0
			elif right >= t[0] and t[1] >= right:
				num_switches += 1
			elif left >= t[0] and t[1] >= left:
				num_switches += 1
	return num_switches

max_window = chrom_length - bene_locus

window = chrom_length
one_sided = []
while window >= max_window:
	window = window/10
	if window <= bene_locus:
		one_sided.append(window)

two_sided = []
while window >= 1e6:
	window = window/10
	two_sided.append(window)

switch_counts = np.zeros(len(one_sided)+1+len(two_sided))

for child in range(8, samples):
	tracts = getTracts(getChildArray(child, ancestor_link))
	left, right = 0, chrom_length

	switch_counts[0] += numSwitches(tracts, left, right)

	for win, i in zip(one_sided, range(len(one_sided))):
		left = bene_locus - win

		switch_counts[i+1] += numSwitches(tracts, left, right)

	for win, i in zip(two_sided, range(len(two_sided))):
		right = bene_locus + win
		left = bene_locus - win

		switch_counts[i+len(one_sided)+1] += numSwitches(tracts, left, right)

switch_counts /= (samples-8)

left = [0.]
right = [chrom_length]

for win in one_sided:
    left.append(bene_locus-win)
    right.append(chrom_length)
for win in two_sided:
    left.append(bene_locus-win)
    right.append(bene_locus+win)

switch_df = pd.DataFrame(np.stack((left, right, switch_counts), axis=-1), columns = ["left", "right", "avgSwitches"])

outname = f"{outfile}_ancestryswitches.csv"
switch_df.to_csv(outname)





