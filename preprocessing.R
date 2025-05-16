if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "IlluminaHumanMethylation450kmanifest",
                       "IlluminaHumanMethylation450kanno.ilmn12.hg19", "DMRcate"))

# ----------------------------------------------------------- #
# Step 1: Download GSE99511 Metadata and Supplementary Files  #
# ----------------------------------------------------------- #
library(GEOquery)
# download manually the files from the GEO repository
# Unzip IDAT files (THIS STATE HAS ONLY PROCESSED ONCE)
untar("GSE99511/GSE99511_RAW.tar", exdir = "C:/Users/putri/Documents/GSE99511/IDATs")

### 
### making sample_name
setwd("C:/Users/putri/Documents/Selected studies/Final selection/GSE99511")
file_list <- list.files("C:/Users/putri/Documents/Selected studies/GSE99511/IDATs", 
                        full.names = TRUE)
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

samplesheet = merge(Sample_Sheet, sample_diseaseStatus, by = "sampleID", all = TRUE)
colnames(samplesheet) = c("Sample_Name", "Sentrix_ID", "Sentrix_Position",
                          "Sample_Group")
write.csv(samplesheet, file="SampleSheet.csv", row.names = FALSE, 
          quote=FALSE)

# ----------------------------------------------------------- #
# Step 2: Read IDAT files and create RGChannelSet             #
# ----------------------------------------------------------- #
library(minfi)

# Locate IDAT files
idat_dir <- "C:/Users/putri/Documents/Selected studies/Final selection/GSE99511"
targets <- read.metharray.sheet(idat_dir)

# Read raw IDAT files
rgSet <- read.metharray.exp(targets = targets)

# Check the sample names and sample sheet
# sampleNames(rgSet)
pData(rgSet) <- targets  # add phenotype data


# ----------------------------------------------------------- #
# Step 3: Preprocess and Normalize                            #
# ----------------------------------------------------------- #
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(tidyr)
library(dplyr)

# Quality control and preprocessing
mSet <- preprocessIllumina(rgSet)  # or use preprocessQuantile / preprocessFunnorm

# Get beta values
Beta <- getBeta(mSet)

# Mapping probes into genes level
# Get the annotation data
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Match CpG IDs in your beta matrix to annotation
common_cpgs <- intersect(rownames(Beta), rownames(anno))
beta_matrix <- Beta[common_cpgs, ]
anno <- anno[common_cpgs, ]

# Create a mapping of CpGs to genes
cpg_gene_map <- strsplit(as.character(anno$UCSC_RefGene_Name), ";")

cpg_to_gene <- data.frame(
  CpG = rep(rownames(beta_matrix), sapply(cpg_gene_map, length)),
  Gene = unlist(cpg_gene_map)
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

ofset = 0.000001
Mvals0 <- log2((gene_level_beta[,-1]+ofset) / (1-gene_level_beta[,-1]+ofset))
rownames(Mvals0) = as.character(gene_level_beta$Gene)

save(gene_level_beta, Mvals0, group, file="dat_genes.RData")

# ----------------------------------------------------------- #
# Step 3: Calculate effect size and variance                  #
# ----------------------------------------------------------- #
library(limma)

load("dat_genes.RData")
load("C:/Users/putri/Documents/Selected studies/Final selection/list_genes.RData")

## take the common probes only
im = match(list_genes$common_genes, rownames(Mvals0))
Mvals = Mvals0[im,]

# #1. VARIANCE
# # Create design matrix
# group0 <- as.factor(pData(mSet)$Sample_Group)  # replace with real group column name
# group = replicate(length(group0), "case")  ## group2 is CIN3+
# group[which(group0=="Normal")] = "control"  ## group1 is Normal
# group = as.factor(group)

design <- model.matrix(~ group)

s2_est <- fit$s2.post #within study variance from limma
n_cont = length(which(group=="control"))
n_case = length(which(group=="case"))
N = n_cont + n_case
J = 1 - (3/(4*N - 9)) #hedges correction

res_table = topTable(fit, number=nrow(Mvals))

## 2. EFFECT SIZE 
teta_ij = xbar_cont_ij = xbar_case_ij = var_cont_ij = var_case_ij = matrix()
for(i in 1:nrow(Mvals)) {
  xbar_cont_ij[i] = mean(as.numeric(Mvals[i,which(group=="control")]), na.rm=TRUE) 
  xbar_case_ij[i] = mean(as.numeric(Mvals[i,which(group=="case")]), na.rm=TRUE)
  
  var_cont_ij[i] = var(as.numeric(Mvals[i,which(group=="control")]), na.rm=TRUE) 
  var_case_ij[i] = var(as.numeric(Mvals[i,which(group=="case")]), na.rm=TRUE) 
  
  xbar = xbar_case_ij[i] - xbar_cont_ij[i]
  teta_ij[i] = (xbar/s2_est[i]) * J
}
summary(teta_ij)

result_MA_GSE99511 = data.frame(teta_ij,s2_est, xbar_cont_ij, xbar_case_ij,
                                var_cont_ij, var_case_ij)
rownames(result_MA_GSE99511) = rownames(Mvals)
head(result_MA_GSE99511)

save(result_MA_GSE99511, file="result_MA_GSE99511.RData")

