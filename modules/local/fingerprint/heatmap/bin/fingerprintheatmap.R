#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)

### Load libraries 
library(VariantAnnotation)
library(dplyr)
library(pheatmap)

# get command arguments
args = commandArgs(trailingOnly = TRUE)
inputVCFs = args[1]
prefix = args[2]

# Read in data
fingerprint_files <- strsplit(inputVCFs, ",")[[1]]
vcfs = lapply(fingerprint_files, function(file) readVcf(file))

# Function to convert the genotyping in solid values
convert_GT <- function( vcf ) {
  x = geno(vcf)$GT
  x[,1] <- c(as.numeric(substr(x,1,1))+as.numeric(substr(x,3,3)))
  return(x)
}

# Visualize the fingerprint data 
df <- as.data.frame(lapply(vcfs, function(x) convert_GT(x)))
df <- df %>% replace(is.na(.), 0)
pdfwidth = 1.5*ncol(df)
pdf(paste0(prefix, "_fingerprintheatmap.pdf"), height = 10, width = pdfwidth)
pheatmap(as.matrix(sapply(df, as.numeric)), fontsize_col = 6, cellwidth = 6, cellheight = 6)
dev.off()





