.libPaths(c(
  "/u/home/v/vturnbil/R/APPTAINER/h2-rstudio_4.5.0",
  .libPaths()
))

library(Mergeomics)

# ---- GET TRAIT FROM SHELL ----
args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0){
  stop("No trait provided!")
}

trait <- args[1]
cat("Running trait:", trait, "\n")

# ---- PATHS ----
base_input  <- "/u/scratch/v/vturnbil/GWAS_post_mdf"
base_output <- "/u/scratch/v/vturnbil/GWAS_MSEA_Mouse"

modfile <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/mouseDE_human.mod"
inffile <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/mouseDE_human.inf"

genfile <- file.path(base_input, trait, "genes.txt")
marfile <- file.path(base_input, trait, "marker.txt")

if(!file.exists(genfile) || !file.exists(marfile)){
  stop("Missing input files for trait: ", trait)
}

outfolder <- file.path(base_output, trait, "MS")
dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

# ---- SETUP MERGEOMICS ----
job.ssea <- list()
job.ssea$label <- trait
job.ssea$folder <- outfolder
job.ssea$genfile <- genfile
job.ssea$locfile <- marfile
job.ssea$marfile <- marfile
job.ssea$modfile <- modfile
job.ssea$inffile <- inffile
job.ssea$permtype <- "gene"
job.ssea$nperm <- 10000
job.ssea$maxoverlap <- 0.33

# ---- RUN ----
job.ssea <- ssea.start(job.ssea)
job.ssea <- ssea.prepare(job.ssea)
job.ssea <- ssea.control(job.ssea)
job.ssea <- ssea.analyze(job.ssea, trim_start=0.005, trim_end=0.995)
job.ssea <- ssea.finish(job.ssea)

cat("Finished:", trait, "\n")