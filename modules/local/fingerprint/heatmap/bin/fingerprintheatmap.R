#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)

### Load libraries 
library(VariantAnnotation)
library(pheatmap)

# get command arguments
args = commandArgs(trailingOnly = TRUE)
inputVCFs = args[1]
prefix = args[2]
#inputVCFs = "~/hpc/pmc_vanboxtel/tools/ToolsVanBox_DEV/ASAP_oatv/demo/output/work/ce/9f9d638fe13b724c26eee180a44f37/PMCAHH1-MSH2KO-MSH2-C25A04SC15C05_dedup.gatk_unifiedgenotyper.vcf.gz, ~/hpc/pmc_vanboxtel/tools/ToolsVanBox_DEV/ASAP_oatv/demo/output/work/b8/c6eb480fa9aff31fbbd5fb48cc7777/PMCAHH1-MSH2KO-MSH2-C25A04_dedup.gatk_unifiedgenotyper.vcf.gz"
#inputDir = "~/hpc/pmc_vanboxtel/projects/_external_projects/GenomeTox/3_Output/FingerPrint/"

# Read in data
fingerprint_files <- strsplit(inputVCFs, ",")[[1]]
print(fingerprint_files)
#test <- readVcf(fingerprint_files[1])




#fingerprint_files <- list.files(path=inputDir, pattern=".vcf$", full.names = T, recursive = T)
#fingerprint_files <- fingerprint_files[1:4]
vcfs = lapply(fingerprint_files, function(file) readVcf(file))

# Function to convert the genotyping in solid values
convert_GT <- function( vcf ) {
  x = geno(vcf)$GT
  x[,1] <- c(as.numeric(substr(x,1,1))+as.numeric(substr(x,3,3)))
  return(x)
}

# Visualize the fingerprint data 
df <- as.data.frame(lapply(vcfs, function(x) convert_GT(x)))
pdf(paste0(prefix, "_fingerprintheatmap.pdf"))
pheatmap(as.matrix(sapply(df, as.numeric)), fontsize_col = 6)
dev.off()
# lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)





