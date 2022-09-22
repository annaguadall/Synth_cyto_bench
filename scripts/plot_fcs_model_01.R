######################################################################################
## PLOT samples in fcs format
## SYNTH_01 project
######################################################################################

## ----Parameters --------------------------------------------------------------------
params <- list(
    log_file = snakemake@log[[1]],
    fcs = snakemake@input[["fcs"]],
    labels = snakemake@input[["labels"]],
    png_56_3 = snakemake@output[["png_56_3"]],
    png_56_19 = snakemake@output[["png_56_19"]],
    png_8_4 = snakemake@output[["png_8_4"]],
    png_3 = snakemake@output[["png_3"]],
    png_56 = snakemake@output[["png_56"]],
    png_19 = snakemake@output[["png_19"]],
    png_8 = snakemake@output[["png_8"]],
    png_4 = snakemake@output[["png_4"]]
)

## ----Save log ----------------------------------------------------------------------
log <- file(params$log_file, open="wt")
sink(log)

## ----Libraries ---------------------------------------------------------------------
library(flowCore)
library(ggcyto)
library(ggridges)

## ----Import data -------------------------------------------------------------------
ff <- read.FCS(params$fcs)
labels <- read.csv(params$labels)[[1]]

## ----Prepare a labelled FlowSet ----------------------------------------------------

# Prepare one FlowFrame per label
for(i in levels(labels)){
    sel <- labels == i
    assign(i, ff[sel,])
}

# Merge all FlowFrames in one FlowSet
fs <- flowSet(B, NK, T4, T8, NKT_NN, NKT_4, NKT_8, U1, U2, U3, U4)
fs@phenoData@data$name <- c("B cells", "NK cells", "T4 cells", "T8 cells", 
                            "CD4-CD8- NKT cells", "CD4+ NKT cells", "CD8+ NKT cells",
                            "U1", "U2", "U3", "U4")

## ----Plot markers expression -------------------------------------------------------
png(params$png_56_3, res=150, width = 1400, height = 1000)
ggcyto(fs, aes(x = "CD56", y = "CD3")) + geom_hex(bins = 50) + xlab("CD56") + ylab("CD3")
dev.off()

png(params$png_56_19, res=150, width = 1400, height = 1000)
ggcyto(fs, aes(x = "CD56", y = "CD19")) + geom_hex(bins = 50) + xlab("CD56") + ylab("CD19")
dev.off()

png(params$png_8_4, res=150, width = 1400, height = 1000)
ggcyto(fs, aes(x = "CD8", y = "CD4")) + geom_hex(bins = 50) + xlab("CD8") + ylab("CD4")
dev.off()

png(params$png_56, res = 150)
ggcyto(fs, aes(x = "CD56", fill = name)) + geom_density() + 
    geom_density_ridges(aes(y = name), alpha = 0.5) + facet_null() + xlab("CD56")
dev.off()

png(params$png_3, res = 150)
ggcyto(fs, aes(x = "CD3", fill = name)) + geom_density() + 
    geom_density_ridges(aes(y = name), alpha = 0.5) + facet_null() + xlab("CD3")
dev.off()

png(params$png_19, res = 150)
ggcyto(fs, aes(x = "CD19", fill = name)) + geom_density() + 
    geom_density_ridges(aes(y = name), alpha = 0.5) + facet_null() + xlab("CD19")
dev.off()

png(params$png_8, res = 150)
ggcyto(fs, aes(x = "CD8", fill = name)) + geom_density() + 
    geom_density_ridges(aes(y = name), alpha = 0.5) + facet_null() + xlab("CD8")
dev.off()

png(params$png_4, res = 150)
ggcyto(fs, aes(x = "CD4", fill = name)) + geom_density() + 
    geom_density_ridges(aes(y = name), alpha = 0.5) + facet_null() + xlab("CD4")
dev.off()

## ----Export data -------------------------------------------------------------------
# Session info 
sessionInfo()

