#! /usr/bin/env Rscript

##############################################################
#Calculates various statistics from tractlenghts.txt
#and ancestryproportions.csv. Specific to outputs from
#CV-admixture_chr1-rates.slim, localancestry_proportions.py,
#and localancestry_tracts.py
#
#Most of this is also specific to the way I've organized
#and named my files. Calculations may be useful, though.
#
#Takes input ancestryproportions.csv file and prints
#the following stats to stdout (so it can be run in
#parallel on cluster. cat output to single file
#after jobs are complete):
#
#median, mean, variance in p1 ancestry tract lenghts,
#percentile rank of beneficial locus ancestry, beneficial
#allele freq, population-wide mean global ancestry,
#95% quantile of local ancestry proportion, 95%
#quantile of tract length, proportion of chromosomes
#with tract length > 75% chromosome 1, proportion of
#chromosomes with tract length > 95% chromosome 1
#
#NOT INCLUDED: proportion of individuals with 
#global ancestry > 0.75. calculated separately in
#globalancestry.py
#
#Requires input csv name similar to following example:
#seed-123_s-0.08_ancestryproportions.csv
#
#and tract lengths file in same directory:
#seed-123_s-0.08_tractlengths.txt
#
#Usage:
#ancestry_stats.R filename.csv
#
##############################################################

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

infile <- commandArgs(trailingOnly=TRUE)[1]

bene_locus <- 159174683.0

#get seed number and selection coeff
pattern <- "seed-(\\d+)_s-(\\d.\\d+)"
pattern_match <- str_match(infile, pattern)
seed <- pattern_match[,2]
s <- pattern_match[,3]

ancestry_prop <- read.csv(infile, header = TRUE)
ancestry_prop %<>% select(-X)
ancestry_prop$P1AncestryProportion <- 1 - ancestry_prop$P2AncestryProportion

#beneficial allele freq
bene_freq <- ancestry_prop[max(which(ancestry_prop$GenomicPosition<=bene_locus)),]$P1AncestryProportion

#mean global ancestry (weighted by interval size)
ancestry_prop$interval <- suppressWarnings(ancestry_prop$GenomicPosition[-1] - ancestry_prop$GenomicPosition)
ancestry_prop <- ancestry_prop[-c(nrow(ancestry_prop)),]
mean_globalancestry <- sum(ancestry_prop$P1AncestryProportion*ancestry_prop$interval/249904549.0)

#95% quantile for ancestry proportion
q95_prop <- quantile(ancestry_prop$P1AncestryProportion, .95, names=FALSE)

#ancestry percentile rank for beneficial locus
bene_rank <- mean(c(ancestry_prop$P1AncestryProportion < bene_freq))

#median, var, mean tract lengths

tract_file <- paste0("./", pattern_match[,1], "_tractlengths.txt")

local_ancestrytracts <- read.table(tract_file, header=TRUE)

median_tractlength <- median(log(local_ancestrytracts$tractLength))

mean_tractlength <- mean(log(local_ancestrytracts$tractLength))

var_tractlength <- var(log(local_ancestrytracts$tractLength))

#95% quantile for tract lengths
q95_tract <-  quantile(log(local_ancestrytracts$tractLength), .95, names=FALSE)

#proportion of individuals with tract length above 75% of the chromosome
length_75_prop <- mean(local_ancestrytracts$tractLength > (249904549.0*.75))
length_95_prop <- mean(local_ancestrytracts$tractLength > (249904549.0*.95))

cat(s, seed, median_tractlength, mean_tractlength, var_tractlength, bene_rank, bene_freq, mean_globalancestry, q95_prop, q95_tract, length_75_prop, paste0(length_95_prop, "\n"), sep="\t")
