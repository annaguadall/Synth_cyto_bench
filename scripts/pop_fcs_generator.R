######################################################################################
## POP GENERATOR
## Normal distribution. Generated samples are exported in fcs format
######################################################################################

## ----Parameters --------------------------------------------------------------------
params <- list(
    log_file = snakemake@log[[1]],
    sample_n = snakemake@input[["sample"]],
    replicate_seed = snakemake@input[["replicate"]],
    neg_mean = snakemake@params[["neg_mean"]],
    neg_sd = snakemake@params[["neg_sd"]],
    pos_mean = snakemake@params[["pos_mean"]],
    pos_sd = snakemake@params[["pos_sd"]],
    phenotypes = snakemake@params[["phenotypes"]],
    labels = snakemake@output[["labels"]],
    fcs = snakemake@output[["fcs"]]
)

## ----Save log ----------------------------------------------------------------------
log <- file(params$log_file, open="wt")
sink(log)

## ----Libraries ---------------------------------------------------------------------
library(flowCore)
library(readxl)

## ----Import data -------------------------------------------------------------------
pheno <- read_xlsx(params$phenotypes)

## ----Prepare parameters ------------------------------------------------------------
# Total number of cells per sample
n <- as.numeric(readLines(params$sample_n)) * 1000

# Seed
seed <- readLines(params$replicate_seed)

# Cell types
cells <- pheno$cells

# Percentage of cells per cell type
pc <- as.numeric(pheno$proportions)

# Phenotype
pheno <- pheno[,-c(1, ncol(pheno))]

# Markers
markers <- colnames(pheno)

# Number of cells per cell type
ns <- c()
for(j in 1:length(pc)){
    ns[j] <- round(n * pc[[j]] / 100)
}

# Final number of cells
final_n <- sum(ns)

## ----POP GENERATION ----------------------------------------------------------------
set.seed(seed)

# CREATING THE DATAFRAME (ITERATE OVER CELL TYPES)
for(l in (1:length(cells))){      # l takes values [1, number of CELL TYPES]
    # For cell type l, repeat its name as many times as the number of cells 
    c <- rep(cells[l], ns[l]) 
    # One column will be created for every marker with intensity values
    for(j in (1:length(markers))){  
        if(pheno[l,j] == "low"){
            p <- rnorm(ns[[l]], mean = params$neg_mean, sd = params$neg_sd)
        }else{
            p <- rnorm(ns[[l]], mean = params$pos_mean, sd = params$pos_sd)
        }
        # all the columns (one per marker) are joined:
        # convert the "c" list as data frame to preserve numeric "p" values after cbind-ing
        c <- cbind(as.data.frame(c), p)  
    }
    if(l==1){
        # Values for the first cell type start to fill the dataframe "pop"
        pop <- c
    }else{
        # Values for the other cell types are joined to "pop"
        pop <- rbind(pop, c)
    }
}

# Name the variables:
names(pop) <- c("cells", markers)

# Randomly sampling
set.seed(42)
x <- sample(1:final_n, final_n, replace = F)

pop <- pop[x,] # Rearranging the dataframe

rownames(pop) <- 1:final_n # Renaming the rows (in order)

# Convert to FlowFrame
ff <- new("flowFrame", exprs = as.matrix(pop[,-1]))

## ----EXPORT -----------------------------------------------------------------------------
# Labels
write.csv(pop$cells, params$labels, row.names = FALSE)

# FlowFrame
write.FCS(ff, params$fcs) 

# Session info 
sessionInfo()
