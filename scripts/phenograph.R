######################################################################################
## PhenoGraph
## 30 nearest neighbours
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
    reduction_time = snakemake@params[["reduction_time"]],
    norm = snakemake@params[["norm"]],
    sample_name = snakemake@params[["sample"]],
    summ = snakemake@output[["summ"]]
)

## ----Save log ----------------------------------------------------------------------
log <- file(params$log_file, open="wt")
sink(log)

## ----Libraries ---------------------------------------------------------------------
library(flowCore)
library(cytofkit)
library(caret)

## ----Functions ---------------------------------------------------------------------
source("scripts/functions.R")

## ----Import data -------------------------------------------------------------------
ff <- read.FCS(params$fcs_file)
labels <- read.csv(params$labels)[[1]]

# Dimensionality reduction time
if(params$reduction_time == FALSE){
    red_user <- NA
    red_elapsed <- NA
}else{
    red_time <- read.csv(params$reduction_time)
    red_user <- red_time[1,2]
    red_elapsed <- red_time[1,3]
}

## ----Clustering --------------------------------------------------------------------
ptm <- proc.time() # measuring CPU time
pred <- Rphenograph(ff@exprs) # default k = 30 nearest neighbours
ptm <- proc.time() - ptm

clustering <- pred$membership

## ----Matching ----------------------------------------------------------------------
matched <- matching(clustering, labels, params$thres) 

## ----Confusion matrix and other performance measurements ---------------------------
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
                   Reduction_user_time = red_user, Reduction_elapsed_time = red_elapsed,
                   Clustering_user_time = ptm[1] + ptm[4], 
                   Clustering_elapsed_time = ptm[3])

write.csv(summ, params$summ)

# Session info 
sessionInfo()



