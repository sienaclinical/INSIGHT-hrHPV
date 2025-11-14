# ============================================================================
# PART 2. Effect Size Calculation for Meta-Analysis
# ============================================================================
# Prerequisites:
# 1. Common genes list: An R object containing genes shared across all studies
#    - Location: list_genes$common_genes
#    - Format: Character vector of gene names
# 
# 2. Methylation data: Log2-transformed beta-values
#    - Format: Matrix or data frame with rownames as gene identifiers
#    - Transformation: log2(beta / (1 - beta))
#    - This is the output from the previous script ("methdat_GSE99511.RData")
# 
# Instructions:
# - Execute this script separately for each study in the integrative analaysis
# - Ensure gene naming conventions are consistent across studies
# 
# Output:
# - Effect size (standardized mean difference, or fold change)
# - Within-study variance
# ============================================================================

library(limma)

load("methdat_GSE99511.RData") ; load("list_genes.RData")

## take the common probes only
im = match(list_genes$common_genes, rownames(Mvals0))
Mvals = Mvals0[im,]

# Create "group"-object (you need an "mset" object, see: "preprocessing_illumina_raw.R")
group0 <- as.factor(pData(mSet)$Sample_Group)  # replace with real group column name
group = replicate(length(group0), "case")  ## group2 is CIN3+
group[which(group0=="Normal")] = "control"  ## group1 is Normal
group = as.factor(group)

# Create design matrix
design <- model.matrix(~ group)

# fit limma model
fit <- lmFit(Mvals, design)
fit <- eBayes(fit)

# WITHIN STUDY VARIANCE (from Limma)
s2_est <- fit$s2.post 

n_cont = length(which(group=="control"))
n_case = length(which(group=="case"))
N = n_cont + n_case

# hedges correction
J = 1 - (3/(4*N - 9)) 

# observe modeling result from limma
res_table = topTable(fit, number=nrow(Mvals))

# EFFECT SIZE 
teta_ij = xbar_cont_ij = xbar_case_ij = var_cont_ij = var_case_ij = matrix()
for(i in 1:nrow(Mvals)) {

  # ---- we need this for latter, to calculate BETWEEN STUDY VARIANCE ------------------ #
  xbar_cont_ij[i] = mean(as.numeric(Mvals[i,which(group=="control")]), na.rm=TRUE) 
  xbar_case_ij[i] = mean(as.numeric(Mvals[i,which(group=="case")]), na.rm=TRUE)
  
  var_cont_ij[i] = var(as.numeric(Mvals[i,which(group=="control")]), na.rm=TRUE) 
  var_case_ij[i] = var(as.numeric(Mvals[i,which(group=="case")]), na.rm=TRUE) 
  # ------------------------------------------------------------------------------------ #
  
  xbar = xbar_case_ij[i] - xbar_cont_ij[i]
  teta_ij[i] = (xbar/s2_est[i]) * J
}
summary(teta_ij)

# save necessary outputs into one R-object
result_GSE99511 = data.frame(teta_ij,s2_est, xbar_cont_ij, xbar_case_ij,
                             var_cont_ij, var_case_ij)
rownames(result_GSE9951) = rownames(Mvals)

save(result_GSE99511, file="result_MA_GSE99511.RData")
