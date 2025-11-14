# ============================================================================
# PART 3: Cross-Study Integration and Meta-Analysis
# ============================================================================
# Purpose:
# Integrate results across multiple studies to perform meta-analysis
#
# Required Inputs:
# 1. Common genes list
#    - Object: list_genes$common_genes
#    - Description: Genes present across all selected studies
# 2. Individual study results
#    - Source: Output from "ind_study.R" script
#    - Required for: All studies included in the meta-analysis
# 3. Study-level statistics
#    - Effect sizes (mean differences per gene)
#    - Within-study variances
#    - Format: One file/object per study
#
# Output:
# - Pooled effect sizes across studies
# - Between-study heterogeneity (I², τ²)
# - Meta-analysis summary statistics
# ============================================================================

library(meta); library(metafor) 

# load datasets
# Directory Structure:
# All outputs are organized by study in separate folders for better management
# 
# Example structure:
# project_root/
# ├── study_GSE99511/
# │   ├── methdat_GSE99511.RData
# ├── study_GSE186835/
# │   ├── methdat_GSE186835.RData
#
# Note: Each study folder contains all intermediate .RData files from each 
# preprocessing and analysis step for reproducibility and troubleshooting.

load("list_genes.RData") #list of common_genes
datasets = c("GSE99511", "GSE143752", "GSE186835", "GSE287994") #GEOID of the selected studies
n_contr = n_case = matrix()
for(i in 1:length(datasets)){
  # load results from "ind_study.R"
  dirX = paste("study_",datasets[i],"/result_", datasets[i],".RData", sep="")
  load(dirX)
  
  dirY = paste(datasets[i],"/methdat_" datasets[i],".RData",sep="")
  load(dirY)
  n_contr[i] = length(which(group=="control"))
  n_case[i] = length(which(group=="case"))
}

# reformatting the data
xbar_contr = cbind(result_GSE99511$xbar_cont_ij, 
                   result_GSE143752$xbar_cont_ij,
                   result_GSE186835$xbar_cont_ij,
                   result_GSE287994$xbar_cont_ij)
xbar_case = cbind(result_GSE99511$xbar_case_ij,
                  result_GSE143752$xbar_case_ij,
                  result_GSE186835$xbar_case_ij,
                  result_GSE287994$xbar_case_ij)
si2_contr = cbind(result_GSE99511$var_cont_ij,
                  result_GSE143752$var_cont_ij,
                  result_GSE186835$var_cont_ij,
                  result_GSE287994$var_cont_ij)
si2_case = cbind(result_GSE99511$var_case_ij,
                 result_GSE143752$var_case_ij,
                 result_GSE186835$var_case_ij,
                 result_GSE287994$var_case_ij)

# ##################################################################### #
# META-ANALYSIS IN 6 DATASETS                                           #
# ##################################################################### #
p = length(list_genes$common_genes)

# Random effect meta-analysis
pval = matrix()
for(g in 1:p){
  meta.model = metacont(n_case,xbar_case[g,],si2_case[g,],
                        n_contr,xbar_contr[g,],si2_contr[g,],
                        method.tau="PM")
  pval[g] = meta.model$pval.random
}
save(pval,file="pval_meta.RData")

# BH corrections
pval.adj = p.adjust(pval,method="BH")
# ##################################################################### #


# ##################################################################### #
# CALCULTE OVERALL EFFECT SIZE                                          #
# ##################################################################### #
teta_ij = cbind(result_GSE143752$teta_ij,
                result_GSE186835$teta_ij,
                result_GSE287994$teta_ij,
                result_GSE99511$teta_ij)
sp2i = cbind(result_GSE143752$s2_est,
             result_GSE186835$s2_est,
             result_GSE287994$s2_est,
             result_GSE99511$s2_est)
rownames(teta_ij) = rownames(sp2i) = rownames(result_GSE143752)

# calculate between-study variance (call the tau2-function first, from "f_tau.R" script)
tau2 = matrix(); 
for(g in 1:nrow(teta_ij)){tau2[g] = estimate.tau2(teta_ij[g,], sp2i[g,], length(datasets), "pm")} 
summary(tau2)

m = matrix(); 
for(g in 1:nrow(teta_ij)){
  y = teta_ij[g,]
  w = 1 / (tau2[g] + sp2i[g,])    
  m[g] = sum(w*y) / sum(w)
}
names(m) = rownames(teta_ij)
results_MA = data.frame(teta = m, tau2, pval.adj)
res = cbind(results_MA, teta_ij, sp2i)

# --------------------------------------------------------- #
# Extracr the differentially methylated genes 
# --------------------------------------------------------- #
alpha = 0.05
idDEG = which(results_MA$pval.adj<alpha)
results_MA_deg = results_MA[idDEG,]
results_MA_deg_sort = results_MA_deg[order(abs(results_MA_deg$teta), decreasing=TRUE),]

idDEG = which(res$pval.adj<alpha)
results_MA_deg = res[idDEG,]
results_MA_deg_sort = results_MA_deg[order(abs(results_MA_deg$teta), decreasing=TRUE),]

