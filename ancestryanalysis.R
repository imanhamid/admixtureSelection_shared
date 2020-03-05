#! /usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

infile <- commandArgs(trailingOnly=TRUE)[1]
coords_file <- "/dscrhome/ih49/admixture_CV/chr1_coords.RData"

chrom_length <- 249904549.0
bene_locus <- 159174683.0

#get seed number and selection coeff
pattern <- "^.+/.+_seed-(\\d+)_s-(\\d.*\\d*)_"
pattern_match <- str_match(infile, pattern)
seed <- pattern_match[,2]
s <- pattern_match[,3]

#ANCESTRY PROPORTION STATS
ancestry_prop <- read.csv(infile, header = TRUE)
ancestry_prop %<>% select(-X)
ancestry_prop$P1AncestryProportion <- 1 - ancestry_prop$P2AncestryProportion

#local ancestry proportion at beneficial locus
bene_prop <- ancestry_prop[max(which(ancestry_prop$GenomicPosition<=bene_locus)),]$P1AncestryProportion

#mean global ancestry (weighted by interval size)
ancestry_prop$interval <- suppressWarnings(ancestry_prop$GenomicPosition[-1] - ancestry_prop$GenomicPosition)
ancestry_prop <- ancestry_prop[-c(nrow(ancestry_prop)),]
mean_globalancestry <- sum(ancestry_prop$P1AncestryProportion*ancestry_prop$interval/chrom_length)

#variance in ancestry proportion
var_prop <- var(ancestry_prop$P1AncestryProportion)

#variance in ancestry proportion in 10 Mb window
tenMb_window_start <- bene_locus - 5000000.0
tenMb_window_end <- bene_locus + 5000000.0

lower <- ancestry_prop[which(ancestry_prop$GenomicPosition>=tenMb_window_start & ancestry_prop$GenomicPosition<=bene_locus),]
upper <- ancestry_prop[which(ancestry_prop$GenomicPosition<=tenMb_window_end & ancestry_prop$GenomicPosition>bene_locus),]

prop_10Mb <- bind_rows(lower, upper)
var_prop_10Mb <- var(prop_10Mb$P1AncestryProportion)

#95% quantile for ancestry proportion
q95_prop <- quantile(ancestry_prop$P1AncestryProportion, .95, names=FALSE)

#ancestry percentile rank for beneficial locus
bene_rank <- mean(c(ancestry_prop$P1AncestryProportion < bene_prop))

#TRACT LENGTH STATS

tract_file <- paste0(pattern_match[,1], "tractlengths.txt")

local_ancestrytracts <- read.table(tract_file, header=TRUE)

#median, var, mean tract lengths
local_ancestrytracts_pos <- local_ancestrytracts[which(local_ancestrytracts$anc=="P1"),]

median_tractlength <- median(log(local_ancestrytracts_pos$tract_length))

mean_tractlength <- mean(log(local_ancestrytracts_pos$tract_length))

var_tractlength <- var(log(local_ancestrytracts_pos$tract_length))

#95% quantile for tract lengths
q95_tract <-  quantile(log(local_ancestrytracts_pos$tract_length), .95, names=FALSE)

#proportion of individuals with tract length above 75% of the chromosome
length_75_prop <- mean(local_ancestrytracts_pos$tract_length > (chrom_length*.75))
length_95_prop <- mean(local_ancestrytracts_pos$tract_length > (chrom_length*.95))

#mean ancestry switches
anc_switches <- mean(local_ancestrytracts$switches)

#proportion of individuals with global >=0.75
local_ancestrytracts$true_global <- c(local_ancestrytracts[which(local_ancestrytracts$anc=="P1"),]$global_ancestry, (1 - local_ancestrytracts[which(local_ancestrytracts$anc=="P2"),]$global_ancestry))

local_ancestrytracts <- local_ancestrytracts[order(local_ancestrytracts$child),]

globalAncestry <- function (child, ancestry_tracts){
  global_ancestry = (ancestry_tracts$true_global[child] + ancestry_tracts$true_global[child+1])/2
  return(global_ancestry)
}

num_samples <- nrow(local_ancestrytracts)
child_indices <- as.data.frame(seq(from=1, to=num_samples, by=2))
globalancestry <- sapply(child_indices, globalAncestry, ancestry_tracts=local_ancestrytracts)
global_75_prop <- sum(c(globalancestry>=0.75))/(num_samples/2)

#unstandardized iHS by ancestry tract
EHH <- function(distance, haplotypes) {
  haplotype_proportions <- sum(c(haplotypes >= distance)) / length(haplotypes)
  EHH <- haplotype_proportions^2
  return(EHH)
}

iHH <- function(distance, EHH){
  suppressMessages(require(pracma))
  trapz_auc <- trapz(distance, EHH)
  return(trapz_auc)
}

load(coords_file)

P1_haplotypes <- c((local_ancestrytracts[which(local_ancestrytracts$anc=="P1"),]$end - bene_locus), (bene_locus- local_ancestrytracts[which(local_ancestrytracts$anc=="P1"),]$start))

P2_haplotypes <- c((local_ancestrytracts[which(local_ancestrytracts$anc=="P2"),]$end - bene_locus), (bene_locus- local_ancestrytracts[which(local_ancestrytracts$anc=="P2"),]$start))

coords$P1_EHH <- sapply(coords$distance,EHH, haplotypes=P1_haplotypes)

coords$P2_EHH <- sapply(coords$distance,EHH, haplotypes=P2_haplotypes)

iHH_P1 <- coords[which(coords$P1_EHH>=0.25),]
iHH_P1 <- iHH_P1[order(iHH_P1$distance),]

iHH_P2 <- coords[which(coords$P2_EHH>=0.25),]
iHH_P2 <- iHH_P2[order(iHH_P2$distance),]

P2_iHH <- iHH(iHH_P2$distance, iHH_P2$P2_EHH)
P1_iHH <- iHH(iHH_P1$distance, iHH_P1$P1_EHH)

iHS_anc <- log(P2_iHH/P1_iHH)

cat(s, seed, bene_rank, bene_prop, mean_globalancestry, var_prop, var_prop_10Mb,
    q95_prop, median_tractlength, mean_tractlength, var_tractlength, q95_tract,
    length_75_prop, length_95_prop, anc_switches, global_75_prop, paste0(iHS_anc, "\n"), sep="\t")