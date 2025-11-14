## IMPORTANT NOTE ##
# Preprocessing Pipeline for Illumina Methylation Array Data
# Purpose: This script demonstrates preprocessing steps for raw IDAT files from Illumina platforms
# Platforms: Illumina HumanMethylation450 (450K) and Illumina MethylationEPIC (850K)
# Dataset: GSE99511 (Illumina HumanMethylation450)

# --------------------- PART 1 ------------------------------------ #
# Preprocess raw datasets downloaded from GEO repository
# You need: the ZIP files (raw) ready on your local repository
# ----------------------------------------------------------------- #

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "IlluminaHumanMethylation450kmanifest",
                       "IlluminaHumanMethylation450kanno.ilmn12.hg19", 
                       "minfi", "DMRcate"))
# DO: set working directory, where the downloaded data is stored
# For example: setwd("C:/Users/Documents/GSE99511")

# ########################################################### #
# Download the dataset and its Supplementary Files            #
# ########################################################### #
library(GEOquery)
# DO: download manually the files from the GEO repository
# Unzip IDAT files 
untar("GSE99511/GSE99511_RAW.tar", exdir = "GSE99511/IDATs")

# -------- making sample_name file -------- #
file_list <- list.files("GSE99511/IDATs", full.names = TRUE)
filenames <- file_list[-c(1:6)]

sample_sheet = data.frame(matrix(NA, nrow = length(filenames), 
                                 ncol = 3))
for(i in 1:length(filenames)){
  filepath <- filenames[i]
  filename <-basename(filepath)

  # Remove the "_Grn.idat" or "_Red.idat" suffix
  core_name <- sub("_[Grn|Red]+\\.idat$", "", filename)
  
  # Split the filename by underscore
  parts <- strsplit(core_name, "_")[[1]]
  
  # Extract parts
  sample_sheet[i,1] <- parts[1]
  sample_sheet[i,2] <- parts[2]
  sample_sheet[i,3] <- parts[3]
}
Sample_Sheet <- sample_sheet[!duplicated(sample_sheet), ] # minus disease class
colnames(Sample_Sheet) = c("sampleID", "Sentrix_ID", "Sentrix_Position")

# get disease class
gset <- getGEO("GSE99511", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
sample_diseaseStatus <- data.frame(disease = gset@phenoData[[1]], 
                                   sampleID = gset@phenoData[[2]])

samplesheet = merge(Sample_Sheet, sample_diseaseStatus, 
                    by = "sampleID", all = TRUE)
colnames(samplesheet) = c("Sample_Name", "Sentrix_ID", "Sentrix_Position",
                          "Sample_Group")
write.csv(samplesheet, file="SampleSheet.csv", row.names = FALSE, 
          quote=FALSE)
# save this csv on the IDATs folder
# you only need to this once

# ########################################################### #
# Read IDAT files and create RGChannelSet                     #
# ########################################################### #
library(minfi)

# Locate IDAT files
idat_dir <- "GSE99511"
targets <- read.metharray.sheet(idat_dir)

# Read raw IDAT files
rgSet <- read.metharray.exp(targets = targets)

# ########################################################### #
# PREPROCESSING AND NORMALIZING                               #
# ########################################################### #
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) #depends on the platform used
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(tidyr)
library(dplyr)

# Quality control and preprocessing
mSet <- preprocessIllumina(rgSet)  # alternatively use preprocessQuantile / preprocessFunnorm

# Get beta values
Beta <- getBeta(mSet)

# ########################################################### #
# Mapping probes into genes level                             #
# ########################################################### #
# Get the annotation data.
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
# Platform-Specific Configuration
# NOTE: Adjust the platform parameter based on your dataset's Illumina technology
# For Illumina MethylationEPIC: use "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"

# Match CpG IDs in your beta matrix to annotation
common_cpgs <- intersect(rownames(Beta), rownames(anno))
beta_matrix <- Beta[common_cpgs, ]
anno <- anno[common_cpgs, ]

# remove missing Regulatory_Feature_Group + unclassified
promoter_probes <- anno[anno$Regulatory_Feature_Group %in% c("Promoter_Associated", 
                                                             "Promoter_Associated_Cell_type_specific",
                                                             "Gene_Associated",
                                                             "Gene_Associated_Cell_type_specific"), ]

# Filter for CpG islands only
cpg_island_probes <- promoter_probes[promoter_probes$Relation_to_Island == "Island", ]

# Match CpG IDs in your beta matrix to annotation
common_cpgs <- intersect(rownames(Beta), rownames(cpg_island_probes))
beta_matrix <- Beta[common_cpgs, ]
anno <- anno[common_cpgs, ]

# 6. Extract gene symbols
cpg_gene_map <- strsplit(as.character(anno$UCSC_RefGene_Name), ";")

# 7. Create mapping dataframe
cpg_to_gene <- data.frame(
  CpG = rep(rownames(beta_matrix), sapply(cpg_gene_map, length)),
  Gene = unlist(cpg_gene_map),
  stringsAsFactors = FALSE
)

# Merge beta values and gene mappings
beta_long <- as.data.frame(beta_matrix)
beta_long$CpG <- rownames(beta_long)
beta_long <- merge(beta_long, cpg_to_gene, by = "CpG")

# Summarize: average beta values by gene per sample
gene_level_beta <- beta_long %>%
  pivot_longer(cols = -c(CpG, Gene), names_to = "Sample", values_to = "Beta") %>%
  group_by(Gene, Sample) %>%
  summarize(Mean_Beta = mean(Beta, na.rm = TRUE)) %>%
  pivot_wider(names_from = Sample, values_from = Mean_Beta)

# log2-transformed beta values 
offset = 0.000001
Mvals0 <- log2((gene_level_beta[,-1]+offset) / (1-gene_level_beta[,-1]+offset))
rownames(Mvals0) = as.character(gene_level_beta$Gene)

save(gene_level_beta, Mvals0, group, file="methdat_GSE99511.RData")
