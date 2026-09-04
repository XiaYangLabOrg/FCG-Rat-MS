#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

.libPaths(c(
  "/u/home/v/vturnbil/R/APPTAINER/h2-rstudio_4.5.0",
  .libPaths()
))

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 2){
  stop("Usage: Rscript run_wKDA.R <network_dir> <deg_rds>")
}

base_dir <- args[1]
deg_file <- args[2]

cat("Base directory:", base_dir, "\n")
cat("DEG file:", deg_file, "\n")

library(biomaRt)
library(dplyr)
library(Mergeomics)

# -------------------------------
# 1. LOAD NETWORK
# -------------------------------

net_path <- file.path(base_dir,"final.network.merged.csv")
net <- read.csv(net_path, stringsAsFactors=FALSE)

colnames(net) <- c("TAIL","HEAD","WEIGHT")

# -------------------------------
# 2. CONVERT TO ENSEMBL IDS
# -------------------------------

mart <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl")

map <- getBM(
  attributes=c("ensembl_gene_id","external_gene_name","rgd_symbol"),
  mart=mart
)

lookup <- c(
  setNames(map$ensembl_gene_id, map$external_gene_name),
  setNames(map$ensembl_gene_id, map$rgd_symbol)
)

convert <- function(x){
  if(grepl("^ENSRNOG",x)) return(x)
  if(x %in% names(lookup)) return(lookup[[x]])
  return(NA)
}

net$TAIL <- sapply(net$TAIL, convert)
net$HEAD <- sapply(net$HEAD, convert)

net <- na.omit(net)
net <- unique(net)
net <- net[net$TAIL != net$HEAD, ]

network_out <- file.path(base_dir,"network_mergeomics.txt")
write.table(net, network_out, sep="\t", row.names=FALSE, quote=FALSE)

# -------------------------------
# 3. LOAD DEG FILE
# -------------------------------

deg <- readRDS(deg_file)

dataset_name <- basename(base_dir)

# ---- FIX: convert folder name to DEG cell type format ----
celltype <- gsub("^MSMC_", "", dataset_name)
celltype <- gsub("_", " ", celltype)
celltype <- sub(" ([0-9]+)$", " (\\1)", celltype)

cat("Detected cell type:", celltype, "\n")

deg_ct <- subset(deg, CellType == celltype)

cat("Number of genes in DEG table:", nrow(deg_ct), "\n")

if(nrow(deg_ct) == 0){
  stop(paste("No DEGs found for cell type:", celltype))
}

# fix Ensembl IDs already in gene column
missing <- is.na(deg_ct$ensembl_gene_id) | deg_ct$ensembl_gene_id==""
is_ens <- grepl("^ENSRNOG", deg_ct$gene)
deg_ct$ensembl_gene_id[missing & is_ens] <- deg_ct$gene[missing & is_ens]

deg_ct <- deg_ct[!is.na(deg_ct$ensembl_gene_id) & deg_ct$ensembl_gene_id!="", ]

# -------------------------------
# 4. MODULE FILE (comparison → DEGs)
# -------------------------------

deg_sig <- deg_ct[deg_ct$adj.P.Val < 0.05, ]

cat("Significant DEGs:", nrow(deg_sig), "\n")

if(nrow(deg_sig) < 10){
  stop("Too few significant DEGs for wKDA modules")
}

module_df <- deg_sig[, c("Comparison","ensembl_gene_id")]
colnames(module_df) <- c("MODULE","NODE")
module_df <- unique(module_df)

modules_out <- file.path(base_dir,"modules.txt")
write.table(module_df, modules_out, sep="\t", row.names=FALSE, quote=FALSE)

# -------------------------------
# 5. GENE SCORE FILE
# -------------------------------

deg_ct$adj.P.Val[deg_ct$adj.P.Val==0] <- 1e-300
deg_ct$SCORE <- -log10(deg_ct$adj.P.Val)

gene_scores <- deg_ct %>%
  group_by(ensembl_gene_id) %>%
  summarise(SCORE=max(SCORE, na.rm=TRUE), .groups="drop")

colnames(gene_scores) <- c("GENE","SCORE")

scores_out <- file.path(base_dir,"gene_scores.txt")
write.table(gene_scores, scores_out, sep="\t", row.names=FALSE, quote=FALSE)

# -------------------------------
# 6. RUN wKDA
# -------------------------------

job.kda <- list()
job.kda$label  <- "wKDA"
job.kda$folder <- base_dir

job.kda$netfile   <- network_out
job.kda$modfile   <- modules_out
job.kda$assocfile <- scores_out

job.kda$edgefactor <- 1
job.kda$depth      <- 1
job.kda$direction  <- 0
job.kda$nperm      <- 2000

moddata <- tool.read(job.kda$modfile)
mod.names <- unique(moddata$MODULE)
moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))),]

subset_file <- file.path(base_dir,"subsetof.supersets.txt")
tool.save(moddata, subset_file)
job.kda$modfile <- subset_file

job.kda <- kda.configure(job.kda)
job.kda <- kda.start(job.kda)
job.kda <- kda.prepare(job.kda)
job.kda <- kda.analyze(job.kda)
job.kda <- kda.finish(job.kda)

job.kda <- kda2cytoscape(job.kda)

cat("wKDA finished successfully\n")