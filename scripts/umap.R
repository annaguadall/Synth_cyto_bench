######################################################################################
## UMAP
## Reduction to 2 dimensions using UMAP
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
library(umap)

## Import data -----------------------------------------------------------------------
ff <- read.FCS(snakemake@input[["fcs"]])

## UMAP dimensionality reduction -----------------------------------------------------
ptm <- proc.time() # measuring CPU time
red <- umap(ff@exprs, random_state = 42)
ptm <- proc.time() - ptm

## EXPORT RESULTS --------------------------------------------------------------------
# UMAP-reduced data saved as FlowFrame
colnames(red$layout) <- c("UMAP_1", "UMAP_2")

ff <- new("flowFrame", exprs = as.matrix(red$layout))
write.FCS(ff, params$red)

# Running time
runtime <- data.frame("Reduction_user_time" = ptm[1] + ptm[4], 
                      "Reduction_elapsed_time" = ptm[3])
rownames(runtime) <- NULL

write.csv(runtime, params$runtime)

# Session info 
sessionInfo()
