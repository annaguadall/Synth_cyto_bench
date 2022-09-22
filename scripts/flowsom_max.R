######################################################################################
## FlowSOM
## Maximum number of clusters = number of expected populations (including ouliers)
######################################################################################

## ----Parameters --------------------------------------------------------------------
params <- list(
    log_file = snakemake@log[[1]],
    fcs_file = snakemake@input[["fcs"]],
    labels = snakemake@input[["labels"]],
    synth = snakemake@params[["synth"]],
    thres = snakemake@params[["thres"]],
    clus = snakemake@output[["clus"]],
    cont = snakemake@output[["cont"]],
    cm = snakemake@output[["cm"]],
    cmm = snakemake@output[["cmm"]],
    f1m = snakemake@output[["f1m"]],
    method = snakemake@params[["method"]],
    norm = snakemake@params[["norm"]],
    sample_name = snakemake@params[["sample"]],
    summ = snakemake@output[["summ"]]
)
## ----Save log ----------------------------------------------------------------------
log <- file(params$log_file, open="wt")
sink(log)

## ----Libraries ---------------------------------------------------------------------
library(flowCore)
library(FlowSOM)
library(caret)

## ----Functions ---------------------------------------------------------------------
source("scripts/functions.R")

## ----Import data -------------------------------------------------------------------
ff <- read.FCS(params$fcs_file)
labels <- read.csv(params$labels)[[1]]

# Number of populations
if(params$synth == "yes") {outliers <- 0}
n_pops <- length(levels(labels)) + outliers

## ----META-CLUSTERING ---------------------------------------------------------------
set.seed(42)

ptm <- proc.time() # measuring CPU time
fs <- FlowSOM(ff, compensate = F,transform = F, scale = F,
              colsToUse = colnames(ff@exprs),  
              nClus = NULL,         # Exact number of clusters for meta-clustering.
              maxMeta = n_pops,     # If NULL, several options will be tried (1:maxMeta)
              seed = 42)
ptm <- proc.time() - ptm

clustering <- fs$metaclustering[fs$FlowSOM$map$mapping[,1]] 

## ----MATCHING ----------------------------------------------------------------------
matched <- matching(clustering, labels, params$thres) 

## ----CONFUSION MATRIX AND OTHER PERFORMANCE MEASUREMENTS ---------------------------
cm <- confusionMatrix(data = matched$clusters[labels != "outliers"],
                      reference = labels[labels != "outliers"])

cm_merged <- confusionMatrix(
    data = matched$merged_clusters[matched$merged_labels != "outliers"],
    reference = matched$merged_labels[matched$merged_labels != "outliers"]
)

mf1 <- mean_f1(cm, cm_merged, matched$merged_labels, matched$count_merged_pops)


## ----Export data -------------------------------------------------------------------
# Clustering results
write.csv(clustering, params$clus, row.names = FALSE)

# Contingency matrix
write.csv(matched$c, params$cont)

# Confusion matrix
write.csv(cm$table, params$cm)

# Merged confusion matrix
write.csv(cm_merged$table, params$cmm)

# F1 matrix
write.csv(matched$f, params$f1m)

# Summary
summ <- data.frame(File = params$sample_name, Norm = params$norm, 
                   n = nrow(ff@exprs), Method = params$method, 
                   Threshold = params$thres, 
                   Clusters = length(table(clustering)), 
                   Partitions = matched$partitions, 
                   Mean_F1 = mf1$mf1, Weighted_mean_F1 = mf1$mf1_w, 
                   Inversed_weights_mean_F1 = mf1$mf1_i_w,
                   Mean_F1_merged = mf1$mf1_m, Corrected_mean_F1_merged = mf1$mf1_m_c,
                   Reduction_user_time = NA, Reduction_elapsed_time = NA,
                   Clustering_user_time = ptm[1] + ptm[4], 
                   Clustering_elapsed_time = ptm[3])

write.csv(summ, params$summ)

# Session info 
sessionInfo()


