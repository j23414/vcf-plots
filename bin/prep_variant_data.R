#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
INFILE <- args[1]

library(tidyverse)
library(readxl)
library(magrittr)
library(Biostrings)

# ==== Testing
INFILE <- "/Users/jchang99/github/j23414/vcf-plots/data/samples_variant_results_longer.xlsx"

CHROM_ORDER <- c("PV062510|PB2", "PV074323|PB1", "PV062508|PA", "PV062513|HA",
                "PV062507|MP", "PV062511|NA", "PV062509|NP", "PV062512|NS")

# Load data
data <- read_excel(INFILE)

# Clean data
complete_data <- data %>%
  select(-IsSynonymous, -AltAminoAcid, -AltCodon, -IsTransition, -AminoAcidChange) %>% # Recompute these based on minor and major variant
  subset(VariantType == "SNP") %>%
  subset(Gene != "PB1-F2") %>%
  mutate(
    Percentage = as.numeric(Percentage),
    CHROM=factor(CHROM, levels=CHROM_ORDER)
  ) %>% # Get all percentages of A,C,G,T at a position
  pivot_wider(
    names_from=ALT,
    values_from = Percentage
  ) %>%
  rowwise() %>%
  mutate( # Sum all percentages in case there are multiple variants at a position
    total=sum(C, A, G, T, na.rm=TRUE)
  )%>%
  ungroup() %>%
  mutate( # Compute the remaining reference percentage based on total
    A=case_when(REF=="A" ~ 1 - total, TRUE ~ A),
    C=case_when(REF=="C" ~ 1 - total, TRUE ~ C),
    G=case_when(REF=="G" ~ 1 - total, TRUE ~ G),
    T=case_when(REF=="T" ~ 1 - total, TRUE ~ T)
  ) %>%
  rowwise() %>%
  mutate( # Specify the major_variant nucleotide
    major_percentage=max(A,C,G,T, na.rm=TRUE),
    consensus=case_when(A==major_percentage ~ "A",
                        G==major_percentage ~ "G",
                        C==major_percentage ~ "C",
                        T==major_percentage ~ "T"),
    SNPCodonPosition = as.numeric(SNPCodonPosition),
    # TODO: Check the other codon positions in case they are in the variant list
    majorCodon =case_when(
      IsGenic ~ paste0(
        substr(RefCodon, 1, SNPCodonPosition),
        consensus,
        substr(RefCodon, SNPCodonPosition + 2, 3)
      ),
      TRUE ~ NA),
    majorAminoAcid = GENETIC_CODE[majorCodon]
  ) %>%
  ungroup()

# Expand back into longer format
expand_data <- complete_data %>%
  pivot_longer(
    cols=c(A,C,G,T),
    names_to = "variant",
    values_to = "minor_percentage"
  ) %>%
  subset(!is.na(minor_percentage)) %>%
  #subset(consensus != variant) %>%
  subset(minor_percentage !=0) %>%
  subset(minor_percentage < 1) %>%
  mutate(
    mutation = paste0(consensus,"->", variant, sep=""),
    minorCodon = case_when(
      IsGenic ~ paste0(
        substr(majorCodon, 1, SNPCodonPosition),
        variant,
        substr(majorCodon, SNPCodonPosition + 2, 3)
      ),
      TRUE ~ NA
    ),
    minorAminoAcid = GENETIC_CODE[minorCodon],
    IsSynonymous = case_when(
      is.na(majorAminoAcid) ~ "Intergenic",
      majorAminoAcid == minorAminoAcid ~ "Synonymous",
      TRUE ~ "Non-Synonymous"
    ), 
    flipped = consensus == REF # diagnostic
  ) %>%
  arrange(Sample)

write_tsv(expand_data, "data/expand_data.tsv")
