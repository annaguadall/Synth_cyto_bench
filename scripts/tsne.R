######################################################################################
## tSNE
## Reduction to 2 dimensions using tSNE
######################################################################################

## ----Parameters --------------------------------------------------------------------
params <- list(
    log_file = snakemake@log[[1]],
    fcs = snakemake@input[["fcs"]],
    red = snakemake@output[["red"]],
    runtime = snakemake@output[["runtime"]]
)

## ----Save log ----------------------------------------------------------------------
log <- file(params$log_file, open="wt")
sink(log)

## Libraries -------------------------------------------------------------------------
library(flowCore)
library(cytofkit)

## Import data -----------------------------------------------------------------------
ff <- read.FCS(snakemake@input[["fcs"]])

## UMAP dimensionality reduction -----------------------------------------------------
ptm <- proc.time() # measuring CPU time
red <-  cytof_dimReduction(data = ff@exprs, method = "tsne", tsneSeed = 42)
ptm <- proc.time() - ptm

## EXPORT RESULTS --------------------------------------------------------------------
# tsne-reduced data saved as FlowFrame
red <- as.data.frame(red)

ff <- new("flowFrame", exprs = as.matrix(red))
write.FCS(ff, params$red) 

# Running time
runtime <- data.frame("Reduction_user_time" = ptm[1] + ptm[4], 
                      "Reduction_elapsed_time" = ptm[3])
rownames(runtime) <- NULL

write.csv(runtime, params$runtime)

# Session info 
sessionInfo()
