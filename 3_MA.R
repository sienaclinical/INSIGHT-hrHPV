# -------------------- PART 3 ---------------------------------------- #
# This script is to integrate information from the selected studies
# You need: list of common_genes, results from "ind_study.R" for all
# your selected studies, mean+variance in individual study
# ------------------------------------------------------------------- #

rm(list=ls()) # clear memory on R - optional

library(meta); library(metafor) 

setwd("C:/Users/Documents/")

# load datasets
load("list_genes.RData") #list of common_genes
datasets = c("GSE99511", "GSE143752", "GSE186835", "GSE287994") #GEOID of the selected studies
n_contr = n_case = matrix()
for(i in 1:length(datasets)){
  # load results from "ind_study.R"
  dirX = paste("C:/Users/Documents/",
  datasets[i],"/result_", datasets[i],".RData", sep="")
  load(dirX)
  
  dirY = paste("C:/Users/Documents/",
               datasets[i],"/dat_genes.RData", sep="")
  load(dirY)
  n_contr[i] = length(which(group=="control"))
  n_case[i] = length(which(group=="case"))
}

# reformatting the data
xbar_contr = cbind(result_GSE99511$xbar_cont_ij, result_GSE143752$xbar_cont_ij,
                      result_GSE186835$xbar_cont_ij, result_GSE287994$xbar_cont_ij)
xbar_case = cbind(result_GSE99511$xbar_case_ij, result_GSE143752$xbar_case_ij,
                      result_GSE186835$xbar_case_ij, result_GSE287994$xbar_case_ij)
si2_contr = cbind(result_GSE99511$var_cont_ij, result_GSE143752$var_cont_ij,
                      result_GSE186835$var_cont_ij, result_GSE287994$var_cont_ij)
si2_case = cbind(result_GSE99511$var_case_ij, result_GSE143752$var_case_ij,
                      result_GSE186835$var_case_ij, result_GSE287994$var_case_ij)

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
teta_ij = cbind(result_GSE143752$teta_ij, result_GSE186835$teta_ij,
                result_GSE287994$teta_ij, result_GSE99511$teta_ij)
sp2i = cbind(result_GSE143752$s2_est, result_GSE186835$s2_est,
             result_GSE287994$s2_est, result_GSE99511$s2_est)
rownames(teta_ij) = rownames(sp2i) = rownames(result_GSE143752)

# calculate between-study variance (call the tau2-function first)
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
head(res)

# --------------------------------------------------------- #
# take the DEG only 
# --------------------------------------------------------- #
alpha = 0.05
idDEG = which(results_MA$pval.adj<alpha)
results_MA_deg = results_MA[idDEG,]
results_MA_deg_sort = results_MA_deg[order(abs(results_MA_deg$teta), decreasing=TRUE),]
head(results_MA_deg_sort,100)

idDEG = which(res$pval.adj<alpha)
results_MA_deg = res[idDEG,]
results_MA_deg_sort = results_MA_deg[order(abs(results_MA_deg$teta), decreasing=TRUE),]
head(results_MA_deg_sort,10)
