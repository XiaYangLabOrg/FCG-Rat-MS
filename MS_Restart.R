# Victoria Turnbill
# 12/17/2025
# FCG Rat MS


################################################################################
# create seurat object using output of cellbender -------------------------
################################################################################

#### Use pak to automatically resolve dependency issues (more intelligent)
# conda install conda-forge::zlib 
# non-zero exit status may be caused by zlib.sh missing which 
# can be found using pak.
# install.packages("pak") 
# pak::pkg_install("igraph")
# install.packages("hdf5r")
#install.packages("scCustomize")

# Install BiocManager if not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# Install scuttle from Bioconductor
#BiocManager::install("scuttle")

library(scuttle)
library(Seurat)
library(pak)
library(igraph)
library(scCustomize)
library(dplyr)
library(ggplot2)
library(stringr)
library(hdf5r)
library(scuttle)

list.files("/u/scratch/v/vturnbil/cellbender_output/MS")

## redo all through MSMC construction to create MSMCXist.RDS from cellbender 
## outputs after adding Xist following Twa et al to check expression of Xist

# Define sample info (folder name + group label)
sample_info <- data.frame(
  sample_id = c('Lane-1__MS1__115-137', 'Lane-1__MS2__119-152', 'Lane-1__MS3__173-235',
               'Lane-1__MS4__149-159', 'Lane-1__MS5__171-129', 'Lane-1__MS6__191-174',
               'Lane-1__MS7__186-197', 'Lane-1__MS8__123-165', 'Lane-1__MS9__190-194'),
  group = c('XYΤ', 'XΧΤ', 'XΧΟ', 
            'XΥΤ', 'XXΤ', 'XXΟ',
            'XΥΤ', 'XXΤ', 'XXΟ'),
  replicate = c('1', '1', '1',
                '2', '2', '2',
                '3', '3', '3'),
  stringsAsFactors = FALSE
)

# Base path where folders are located
base_path <- "/u/scratch/v/vturnbil/cellbender_output/MS"

# Create a list to hold individual Seurat objects
seurat_list <- list()

# Loop through each sample
for (i in 1:nrow(sample_info)) {
  sid <- sample_info$sample_id[i]
  group <- sample_info$group[i]
  replicate <- sample_info$replicate[i]
  
  # Path to CellBender output (assumes 10x-style .h5)
  h5_path <- file.path(base_path, sid, paste0(sid,".cellbender_filtered.h5"))
  
  # Read and create Seurat object
  counts <- Read_CellBender_h5_Mat(h5_path)
  seu <- CreateSeuratObject(counts = counts, project = sid, min.cells = 3) #, min.features = 200
  
  # Add metadata
  seu$sample <- sid
  seu$group <- group
  seu$replicate <- replicate
  
  # Add metadata
  seu$sample <- sid
  seu$group <- group
  seu$replicate <- replicate
  
  # Store
  seurat_list[[sid]] <- seu
}

# Merge all samples into one Seurat object
combined <- merge(seurat_list[[1]], y = seurat_list[-1], 
                  add.cell.ids = names(seurat_list), project = "MS")

# Add % mito content
combined[["percent.mt"]] <- 
  PercentageFeatureSet(combined, pattern = "^Mt-")  # or ^mt- for mouse

# save seurat object as HP1 ----------------------------------------------------
# saveRDS(combined,"/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS1.RDS")
# combined <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS1.RDS")
# ------------------------------------------------------------------------------


################################################################################
# QC Filtering of Seurat Object
################################################################################


table(combined$orig.ident)

Layers(combined[["RNA"]])# before joinlayers

combined <- JoinLayers(combined)

Layers(combined[["RNA"]])# after joinlayers

1:9

paste("^Gm",1:9)

# if working with human data, change the appropriate prefixes
# Add QC metrics to Seurat object
seuratObject <- Add_Cell_QC_Metrics(
  object = combined,
  species = "rat",   # "rat", "mouse", or "human"
  feature_patterns = list(
    percent.mito = "^Mt-",                    # Mitochondrial genes
    percent.ribo = c("^Rps", "^Rpl"),        # Ribosomal genes
    percent.pred = c("Gm1","Gm2"),           # Predicted genes
    percent.Hb = c("^Hba-","^Hbb-","^Hbq")   # Hemoglobin genes
  )
)

seuratObject[["percent.mito"]] <- PercentageFeatureSet(
  seuratObject,
  pattern = "^Mt-"
)

seuratObject[["percent.ribo"]] <- PercentageFeatureSet(
  seuratObject,
  pattern = "^Rps|^Rpl"
)

seuratObject[["percent.pred"]] <- PercentageFeatureSet(
  seuratObject,
  features = grep("^Gm", rownames(seuratObject), value = TRUE)
)

seuratObject[["percent.Hb"]] <- PercentageFeatureSet(
  seuratObject,
  pattern = "^Hba-|^Hbb-|^Hbq"
)

p <- VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/violin.pdf",
  plot = p,
  dpi = 300
)

library(patchwork)

# Violin plots for QC metrics
p_list <- VlnPlot(
  seuratObject,
  features = c("nFeature_RNA","nCount_RNA","percent.mito","percent.ribo","percent.pred"),
  pt.size = 0,
  combine = FALSE
)

# Add font size to each plot
p_list <- lapply(p_list, function(x) {
  x + theme(legend.position = "none")
})

# Combine
p <- wrap_plots(p_list, ncol = 3)

# Save
ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/violin.pdf",
  plot = p,
  width = 20,
  height = 10,
  dpi = 300
)

# Scatter plots for QC metrics
p <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/features.pdf",
  plot = p,
  dpi = 300
)

quantile(seuratObject$nFeature_RNA, probs = c(0.02, 0.95))
quantile(seuratObject$nFeature_RNA, probs = c(0.1, 0.95))
quantile(seuratObject$nFeature_RNA, probs = c(0.2, 0.99))

quantile(seuratObject$nCount_RNA, probs = c(0.02, 0.95))
quantile(seuratObject$nCount_RNA, probs = c(0.15, 0.9))
quantile(seuratObject$nCount_RNA, probs = c(0.2, 0.92))

filtered_sample <- subset(
  x = seuratObject,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3500 &
    seuratObject$percent_mito < 0.01 &
    #     seuratObject$percent_ribo < 0.05 &
    #     seuratObject$percent.pred < 0.01 &
    nCount_RNA > 200 &
    nCount_RNA < 4500 ## set according to sequencing depth. 20k per cell
)

p <- VlnPlot(filtered_sample, features = c("nFeature_RNA","nCount_RNA","percent.mito","percent.ribo","percent.pred"), 
             ncol = 3, pt.size = 0)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/violin_filtered.pdf",
  plot = p,
  width = 4*5,
  dpi = 300
)

p <- FeatureScatter(filtered_sample, feature1 = "nCount_RNA", 
                    feature2 = "nFeature_RNA")

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/features_filtered.pdf",
  plot = p,
  dpi = 300
)

table(filtered_sample$orig.ident)

num_cells <- nrow(filtered_sample@meta.data)
print(num_cells)

# save as MS2 ------------------------------------------------------------------
# saveRDS(filtered_sample, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS2.RDS")
# filtered_sample <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS2.RDS")
# ------------------------------------------------------------------------------


################################################################################
# Doublet Removal (DoubletFinder)
################################################################################

# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

library(DoubletFinder)
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(tibble)

options(bitmapType = "cairo")
pdf(NULL)

# run_doubletfinder_custom
# Here's the official rule of thumb from 10x: 
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled

## add more rows of recovered numbers
original.rate <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076, 0.1, 0.16,0.2,0.24,0.3),
                            'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000, 25600, 32000,40000,48000,60000),
                            'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 16000,20000,25000,30000,38000))


# create function
# run_doubletfinder_custom runs Doublet_Finder() and returns a dataframe with the cell IDs and a column with either 'Singlet' or 'Doublet'
run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL){
  # for debug
  #seu_sample_subset <- samp_split[[1]]
  # Print sample number
  print(paste0("Sample ", unique(seu_sample_subset[['orig.ident']]), '...........')) ## change to "orig.ident"
  
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076, 0.1, 0.16,0.2,0.24,0.3),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000, 25600, 32000,40000,48000,60000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 16000,20000,25000,30000,38000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # Pre-process seurat object with standard seurat workflow --- 
  sample <- NormalizeData(seu_sample_subset)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  
  # Finish pre-processing with min_pc
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
  sample <- FindClusters(object = sample, resolution = 0.1)
  
  # pK identification (no ground-truth) 
  #introduces artificial doublets in varying props, merges with real data set and 
  # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
  # provides a list of the proportion of artificial nearest neighbours for varying
  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  ## Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # run DoubletFinder
  sample <- doubletFinder(seu = sample, #reuse.pANN = FALSE,## add reuse.pANN = TRUE 2025-05-09
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  # change name of metadata column with Singlet/Doublet information
  # new columns were added to seurat metadata. pANN_0.25_0.005_3981	DF.classifications_0.25_0.005_3981
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
  colnames(sample@meta.data)[grepl('pANN.*', colnames(sample@meta.data))] <- "pANN_score"
  
  
  # Subset and save
  # head(sample@meta.data['doublet_finder'])
  # singlets <- subset(sample, doublet_finder == "Singlet") # extract only singlets
  # singlets$ident
  double_finder_res <- sample@meta.data[c('doublet_finder','pANN_score')] # get the metadata column with singlet, doublet info
  double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
  return(double_finder_res)
}


# DoubletFinder should be run on a per sample basis

library(tidyverse)


samp_split <- SplitObject(filtered_sample, split.by = "orig.ident")

# Get Doublet/Singlet IDs by DoubletFinder()

samp_split1 <- lapply(samp_split, run_doubletfinder_custom) 

dev.off()

# get singlet/doublet assigned to each of the cell IDs (each element of the 
# list is a different sample)


# You can also set the multiplet rate manually:
# samp_split <- # samp_split <-lapply(samp_split, run_doubletfinder_custom, multiplet_rate = 0.01)


sglt_dblt_metadata <- data.frame(bind_rows(samp_split1)) # merge to a single dataframe
rownames(sglt_dblt_metadata) <- sglt_dblt_metadata$row_names # assign cell IDs to row names to ensure match
sglt_dblt_metadata$row_names <- NULL
head(sglt_dblt_metadata)
filtered_sample.doulbeFinder <- AddMetaData(filtered_sample, sglt_dblt_metadata, col.name = c('doublet_finder','pANN_score'))


table(filtered_sample.doulbeFinder$sample,filtered_sample.doulbeFinder$doublet_finder)


# save as MS3.RDS --------------------------------------------------------------
# saveRDS(filtered_sample.doulbeFinder,"/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS3.RDS")
# filtered_sample.doulbeFinder <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS3.RDS")
# ------------------------------------------------------------------------------


# Check how doublets singlets differ in QC measures per sample.
p <- VlnPlot(filtered_sample.doulbeFinder, group.by = 'orig.ident', split.by = "doublet_finder",
             features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo", "percent.Hb"), 
             ncol = 2, pt.size = 0) + theme(legend.position = 'right')

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/violin_doublet.pdf",
  plot = p,
  height = 4 *3,
  width = 5 *3,
  dpi = 300
)

# Removing doublets
# If all is ok subset and remove doublets
MS.single <- subset(filtered_sample.doulbeFinder, doublet_finder == 'Singlet')

num_cells <- nrow(MS.single@meta.data)
print(num_cells)

# save as MS4 ------------------------------------------------------------------
# saveRDS(MS.single, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS4.RDS")
# MS.single <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS4.RDS")
# ------------------------------------------------------------------------------

## check the distribution of doublets

# use the TransferData() from seurat for cell annotation
ref <- NormalizeData(filtered_sample.doulbeFinder)
ref <- FindVariableFeatures(ref)
ref <- ScaleData(ref)
ref <- RunPCA(ref) #default 50 PC

ref <- FindNeighbors(ref, reduction = "pca", dims = 1:30)
ref <- FindClusters(ref)
MS.ref <- RunUMAP(ref, dims = 1:30)

Idents(MS.ref) <- "doublet_finder"
p <- DimPlot(MS.ref,reduction = "umap",  label = FALSE)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/umap_doublet.pdf",
  plot = p,
  dpi = 300
)

p <- FeaturePlot(MS.ref, features = "nCount_RNA")

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/features_doublet.pdf",
  plot = p,
  dpi = 300
)


################################################################################
# Doublet Removal via Scrublet
################################################################################

# remotes::install_github("Moonerss/scrubletR")
# BiocManager::install("glmGamPoi")

# install.packages("R.utils", lib = "/u/home/v/vturnbil/R/APPTAINER/h2-rstudio_4.5.0")
# install.packages("R.methodsS3", lib = "/u/home/v/vturnbil/R/APPTAINER/h2-rstudio_4.5.0")
# install.packages("R.oo", lib = "/u/home/v/vturnbil/R/APPTAINER/h2-rstudio_4.5.0")

# install.packages(c("usethis", "gitcreds"))
# usethis::create_github_token()
# gitcreds::gitcreds_set()

# remotes::install_github("satijalab/seurat-wrappers",
#                         lib = "/u/home/v/vturnbil/R/APPTAINER/h2-rstudio_4.5.0",
#                         dependencies = TRUE,
#                         force = TRUE)


library(glmGamPoi) 
library(Seurat)
library(SeuratWrappers)
library(scrubletR)
library(reticulate)


seu <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS3.RDS")
# with pANN_score


# run in linux terminal:
# conda create -n py310 python=310
# conda activate py310
# conda install -c conda-forge numpy pandas scanpy

# conda activate x2
# conda install -c conda-forge python-annoy -y
# pip install --no-cache-dir scrublet==0.2.3
# conda install -c conda-forge numpy -y


python_path <- "/u/project/xyang123/vturnbil/packages/miniconda3/envs/x2/bin"


# Explicitly set the Python path
reticulate::use_python(python_path, required = TRUE)
use_python(python_path, required = TRUE)

# Check Python configuration
reticulate::py_config()


samples <- unique(seu$orig.ident)
seu$scrublet_doublet <- NA
seu$scrublet_score <- NA


for (sample_id in samples) {
  print(paste("Processing sample:", sample_id))
  # Running scrublet
  scrublet_out <- scrublet_R(seu[, seu$orig.ident == sample_id],python_home = python_path)
  
  # Store the scrublet results back into the Seurat object's metadata.
  cells_in_sample <- colnames(seu)[seu$orig.ident == sample_id]
  seu@meta.data[cells_in_sample, "scrublet_doublet"] <- as.factor(scrublet_out$predicted_doublets)
  seu@meta.data[cells_in_sample, "scrublet_score"] <- scrublet_out$doublet_scores
}


# save as MS5.RDS --------------------------------------------------------------
# saveRDS(seu,"/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS5.RDS") # with pANN_score
# filtered_sample.doulbeFinder_Scrublet <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS5.RDS")
# ------------------------------------------------------------------------------

table(seu$scrublet_doublet)
table(seu$doublet_finder)

library(ggplot2)

options(repr.plot.width = 15, repr.plot.height = 10)
p <- ggplot(seu@meta.data, aes(x = scrublet_score, fill = orig.ident)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  facet_wrap(~ orig.ident, scales = "free_y") +
  theme_bw() +
  labs(title = "Scrublet Doublet Scores per Sample")

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/scrublet.pdf",
  plot = p,
  dpi = 300
)

## check the distribution of doublets

# use the TransferData() from seurat for cell annotation
ref <- NormalizeData(seu)
ref <- FindVariableFeatures(ref)
ref <- ScaleData(ref)
ref <- RunPCA(ref) #default 50 PC

ref <- FindNeighbors(ref, reduction = "pca", dims = 1:30)
ref <- FindClusters(ref)
MS.ref <- RunUMAP(ref, dims = 1:30)

Idents(MS.ref) <- "scrublet_doublet"
p <- DimPlot(MS.ref,reduction = "umap",  label = FALSE)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/umap_scrublet.pdf",
  plot = p,
  dpi = 300
)

p <- FeaturePlot(MS.ref, features = "nCount_RNA")

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/features_scrublet.pdf",
  plot = p,
  dpi = 300
)


library(dplyr)

# find appropriate threshold values for each sample:
range(filter((seu[, seu$orig.ident == "Lane-1__MS1__115-137"][[]]),scrublet_doublet=="2")$scrublet_score)

thresholds <- c("Lane-1__MS1__115-137" = 0.52,
                "Lane-1__MS2__119-152" = 0.76,
                "Lane-1__MS3__173-235" = 0.36,
                "Lane-1__MS4__149-159" = 0.74,
                "Lane-1__MS5__171-129" = 0.72,
                "Lane-1__MS6__191-174" = 0.51,
                "Lane-1__MS7__186-197" = 0.75,
                "Lane-1__MS8__123-165" = 0.42,
                "Lane-1__MS9__190-194" = 0.4)

unique(seu$orig.ident)
filtered_scrublet <- subset(seu,
                            subset = !((orig.ident == "Lane-1__MS1__115-137" & scrublet_score > thresholds["Lane-1__MS1__115-137"]) |
                                         (orig.ident == "Lane-1__MS2__119-152" & scrublet_score > thresholds["Lane-1__MS2__119-152"]) |
                                         (orig.ident == "Lane-1__MS3__173-235" & scrublet_score > thresholds["Lane-1__MS3__173-235"]) |
                                         (orig.ident == "Lane-1__MS4__149-159" & scrublet_score > thresholds["Lane-1__MS4__149-159"]) |
                                         (orig.ident == "Lane-1__MS5__171-129" & scrublet_score > thresholds["Lane-1__MS5__171-129"]) |
                                         (orig.ident == "Lane-1__MS6__191-174" & scrublet_score > thresholds["Lane-1__MS6__191-174"]) |
                                         (orig.ident == "Lane-1__MS7__186-197" & scrublet_score > thresholds["Lane-1__MS7__186-197"]) |
                                         (orig.ident == "Lane-1__MS8__123-165" & scrublet_score > thresholds["Lane-1__MS8__123-165"]) |
                                         (orig.ident == "Lane-1__MS9__190-194" & scrublet_score > thresholds["Lane-1__MS9__190-194"]))
)


table(filtered_scrublet$sample,filtered_scrublet$scrublet_doublet)

#aggressive--->less singlet
filtered_scrublet$union_doublet <- ifelse(filtered_scrublet$doublet_finder == "Doublet" | filtered_scrublet$scrublet_doublet == "2", "Doublet", "Singlet")

table(filtered_scrublet$sample,filtered_scrublet$union_doublet)

#conservative
filtered_scrublet$combined_doublet <- ifelse(filtered_scrublet$doublet_finder == "Doublet" & filtered_scrublet$scrublet_doublet == "2", "Doublet", "Singlet")

table(filtered_scrublet$sample,filtered_scrublet$combined_doublet)

# used aggressive version ^^^

num_cells <- nrow(filtered_scrublet@meta.data)
print(num_cells)

# Save as MS6.RDS --------------------------------------------------------------
# saveRDS(filtered_scrublet,"/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS6.RDS")
# filtered_sample.doulbeFinder_Scrublet_union <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS6.RDS")
# ------------------------------------------------------------------------------

################################################################################
# create new seurat object with only singlets (as defined by both DoubletFinder
# and Scrublet)
################################################################################

MS <- subset(filtered_sample.doulbeFinder_Scrublet_union, subset=union_doublet=="Singlet")

num_cells <- nrow(MS@meta.data)
print(num_cells)

################################################################################
# Seurat QC and Integration
################################################################################

library(Seurat)
library(SeuratObject) 
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(harmony)
library(SeuratWrappers)


MS <- NormalizeData(MS)
MS <- FindVariableFeatures(MS, selection.method = "vst", nfeatures = 2000)
MS <- ScaleData(MS)
MS <- RunPCA(MS) #default 50 PC

p <- ElbowPlot(MS, ndims = 50)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/elbow.pdf",
  plot = p,
  dpi = 300
)


# using 40 PCs based on elbow plot

MS <- FindNeighbors(MS, reduction = "pca", dims = 1:40)
MS <- FindClusters(MS)
MS <- RunUMAP(MS, dims = 1:40)


Idents(MS) <- "RNA_snn_res.0.8"
p <- DimPlot(MS,reduction = "umap",  label = TRUE)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/umap.pdf",
  plot = p,
  dpi = 300
)

# discovered typos when making groups originally (Greek vs English alphabet). Fix here:
unique(MS$group)

MS$group <- MS$group %>%
  stringr::str_replace_all(c(
    "Χ" = "X",  # Greek chi → X
    "Τ" = "T",  # Greek tau → T
    "Ο" = "O",  # Greek omicron → O
    "Υ" = "Y"   # Greek upsilon → Y
  ))

unique(MS$group)
# should now give only 3 groups, e.g. "XYT", "XXT", "XXO"



# visualize according to group prior to integration

Idents(MS) <- "group"
p <- DimPlot(MS,reduction = "umap",  label = FALSE)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/umap_group.pdf",
  plot = p,
  dpi = 300
)


# assign batch information from GSU and GZ sample handling
library(stringr)

replace_vec <- c('Lane-1__MS1__115-137'='batch1', 'Lane-1__MS2__119-152'='batch1', 'Lane-1__MS3__173-235'='batch1',
                 'Lane-1__MS4__149-159'='batch2', 'Lane-1__MS5__171-129'='batch2', 'Lane-1__MS6__191-174'='batch2',
                 'Lane-1__MS7__186-197'='batch3', 'Lane-1__MS8__123-165'='batch3', 'Lane-1__MS9__190-194'='batch3')

MS$batch <- str_replace_all(MS$sample, replace_vec)

unique(MS$batch)

# visualize according to batch prior to integration
Idents(MS) <- "batch"
p <- DimPlot(MS,reduction = "umap",  label = FALSE)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/umap_batch.pdf",
  plot = p,
  dpi = 300
)

# visualize prior to integration --> need integration
p <- FeaturePlot(MS, features = "nCount_RNA")

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/umap_group.pdf",
  plot = p,
  dpi = 300
)

# perform integration with Harmony

MS <- RunHarmony(MS, "batch")

# Harmony converged after 3 iterations


# downstream
MS <- RunUMAP(MS, reduction = "harmony", dims = 1:40)
MS <- FindNeighbors(MS, reduction = "harmony", dims = 1:40)
MS <- FindClusters(MS)
MS <- JoinLayers(MS)


# visualize after integration
Idents(MS) <- "RNA_snn_res.0.8"
p <- DimPlot(MS,reduction = "umap",  label = TRUE)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/umap_int.pdf",
  plot = p,
  dpi = 300
)

# visualize by group after integration
Idents(MS) <- "group"
p <- DimPlot(MS,reduction = "umap",  label = FALSE)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/umap_group_int.pdf",
  plot = p,
  dpi = 300
)

# visualize by batch after integration
Idents(MS) <- "batch"
p <- DimPlot(MS,reduction = "umap",  label = FALSE)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/umap_batch_int.pdf",
  plot = p,
  dpi = 300
)

# Save as MS7.RDS --------------------------------------------------------------
# saveRDS(MS,"/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS7.RDS")
# MS <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS7.RDS")
# ------------------------------------------------------------------------------

################################################################################
# verify genotypes
################################################################################

counts_ut <- GetAssayData(MS, assay = "RNA", slot = "counts")["Uty", ]

tapply(counts_ut, MS$sample, sum)

# 1, 4, 7

counts_ut <- GetAssayData(MS, assay = "RNA", slot = "counts")["ENSRNOG00000065796", ]

tapply(counts_ut, MS$sample, sum)

# 2, 3, 5, 6, 8, 9

# do with Xist after adding Xist following Twa et al

################################################################################
# Begin Cell Type Annotation
################################################################################


# convert rat genes into mouse orthologs in preparation for MapMyCells using GeneOrthology

MSmouse <- MS


# download ortholog table: https://github.com/AllenInstitute/GeneOrthology/blob/main/csv/mammalian_orthologs_20231113.csv


# Read in a conversion table generated from GeneOrthology
convert <- read_csv("/u/project/xyang123/xyang123-NOBACKUP/turnbill/Tools/MapMyCells/mammalian_orthologs_20231113.csv")

# Install (if needed), then load R library with query data set
# --If you are using your own, read it in here to a variable called "dataIn" with genes as row names
if(!is.element("hodge2019data",.packages(all.available = TRUE)))  # Install data package if not already installed
  devtools::install_github("AllenInstitute/hodge2019data")
library(hodge2019data)


mat_norm <- MSmouse[["RNA"]]$data  # sparse log-normalized matrix
dataIn <- mat_norm

colnames(convert)

# convert rat genes from a mix of EnsemblIDs and gene symbols to just EnsbmblIDs

symbol_map <- convert[
  !is.na(convert$Rat_Symbol) & !is.na(convert$Rat_EnsemblID),
  c("Rat_Symbol", "Rat_EnsemblID")
]

new_names <- rownames(dataIn)

symbol_idx <- new_names %in% symbol_map$Rat_Symbol

new_names[symbol_idx] <- symbol_map$Rat_EnsemblID[
  match(new_names[symbol_idx], symbol_map$Rat_Symbol)
]

rownames(dataIn) <- new_names

# check for duplicated genes (two mouse genes for one rat gene)
any(duplicated(symbol_map$Rat_Symbol))

# Do the conversion by EnsemblID
convert_by_EnsemblID <- convert[!(is.na(convert$Rat_EnsemblID)|is.na(convert$Mouse_EnsemblID)),c("Rat_EnsemblID","Mouse_EnsemblID")] # Remove NAs from conversion table
convert_by_EnsemblID <- convert_by_EnsemblID[is.element(convert_by_EnsemblID$Rat_EnsemblID,rownames(dataIn)),] # Remove genes not in data matrix


dataOut <- dataIn[match(convert_by_EnsemblID$Rat_EnsemblID,rownames(dataIn)),] # Subset data to include only genes with mouse othologs
rownames(dataOut) <- convert_by_EnsemblID$Mouse_EnsemblID # Convert gene names to mouse

# Output the new data matrix

# these two steps will take forever

write.csv(dataOut, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/MapMyCells/converted_data_matrixXist.csv")
# Note: for output in a format compatible with MapMyCells, see documentation here:
#  https://portal.brain-map.org/explore/file-requirements-and-limits

converted_data_matrix <- read_csv("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/MapMyCells/converted_data_matrixXist.csv")
# -----------------------------------------------------------------------------

# convert csv file into new MSmouse seurat object and save

library(Seurat)
library(Matrix)

# Step 2: Convert to matrix with rownames
mat <- as.matrix(converted_data_matrix[,-1])   # numeric part
rownames(mat) <- converted_data_matrix[[1]]    # first column = gene names

# Step 3: Convert to sparse matrix to save memory
mat_sparse <- Matrix(mat, sparse = TRUE)

# Step 4: Create a new Seurat object
MSmouse <- CreateSeuratObject(counts = mat_sparse, assay = "RNA")

# Step 5: Check the object
MSmouse
nrow(MSmouse)  # number of genes
ncol(MSmouse)  # number of cells

# save as seurat object with mouse genes ---------------------------------------
# saveRDS(MSmouse, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/MapMyCells/MSmouse.RDS")
# MSmouse <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/MapMyCells/MSmouse.RDS")
# -----------------------------------------------------------------------------

# convert MSmouse to h5ad file for use with MapMyCells

## install myfork, which compatible with latest anndata 
# remotes::install_github("zqfang/MuDataSeurat", force = T)
library(MuDataSeurat)

# (optional) step 1: Slim down a Seurat object. So you get raw counts, lognorm counts
# Make a copy of the RNA assay
rna_assay <- MSmouse[["RNA"]]

# Remove scale.data if it exists
if ("scale.data" %in% names(rna_assay@layers)) {
  rna_assay@layers[["scale.data"]] <- NULL
}

# Put the assay back
MSmouse[["RNA"]] <- rna_assay

# Now DietSeurat
seu <- DietSeurat(
  object = MSmouse,
  assays = "RNA",
  features = rownames(MSmouse),
  layers = c("counts", "data"),
  graphs = c("RNA_nn", "RNA_snn"),
  dimreducs = c("pca", "umap"),
  misc = TRUE
)

# for seurat v5, need to JoinLayer first
DefaultAssay(MSmouse) = "RNA"
MSmouse <- JoinLayers(MSmouse)
# single modality
MuDataSeurat::WriteH5AD(MSmouse, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/MapMyCells/MSmouse.h5ad", assay="RNA")
# multi modality, ATAC+RNA, CITE-seq 
# MuDataSeurat::WriteH5MU(MS, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/MapMyCells/MSmouse.h5mu")

#### to load MSmouseH5AD file for use with MapMyCells
# MSmouseH5AD <- MuDataSeurat::ReadH5AD("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/MapMyCells/MSmouse.h5ad")
# -----------------------------------------------------------------------------

# ran successfully in MapMyCells

# load .csv output from MapMyCells for MSmouseH5AD

mapmycells <- read.csv("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/MapMyCells/MSmouse_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1767828657113/MSmouse_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1767828657113.csv",
                       skip = 4,               # skip first 4 metadata lines
                       stringsAsFactors = FALSE,
                       row.names = 1           # first column 'cell_id' becomes rownames
)


# Example: add class-level annotation
MSmouse$MapMyCells_Class <- mapmycells[colnames(MSmouse), "class_name"]

# Example: add subclass-level annotation
MSmouse$MapMyCells_Subclass <- mapmycells[colnames(MSmouse), "subclass_name"]

head(MSmouse@meta.data)
table(MSmouse$MapMyCells_Class)
table(MSmouse$MapMyCells_Subclass)


# save MSmouse2 after adding metadata from MapMyCells --------------------------
# saveRDS(MSmouse, file = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/MapMyCells/MSmouse2.rds")
# MSmouse <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/MapMyCells/MSmouse2.rds")
# ------------------------------------------------------------------------------

# add mouse gene names into metadata in MS seurat object


library(Seurat)

# Read in a conversion table generated from GeneOrthology
convert <- read_csv("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/HP/MapMyCells/mammalian_orthologs_20231113.csv")

library(Seurat)

# convert rat genes from a mix of EnsemblIDs and gene symbols to just EnsbmblIDs

mat_norm <- MS[["RNA"]]$data  # sparse log-normalized matrix
dataIn <- mat_norm

colnames(convert)

symbol_map <- convert[
  !is.na(convert$Rat_Symbol) & !is.na(convert$Rat_EnsemblID),
  c("Rat_Symbol", "Rat_EnsemblID")
]

Rat_EnsemblID <- rownames(dataIn)

symbol_idx <- Rat_EnsemblID %in% symbol_map$Rat_Symbol

Rat_EnsemblID[symbol_idx] <- symbol_map$Rat_EnsemblID[
  match(Rat_EnsemblID[symbol_idx], symbol_map$Rat_Symbol)
]

rownames(dataIn) <- Rat_EnsemblID

# check for duplicated genes (two mouse genes for one rat gene)
any(duplicated(symbol_map$Rat_Symbol))


# Step 2: Prepare mapping
convert_by_EnsemblID <- convert[!(is.na(convert$Rat_EnsemblID) | is.na(convert$Mouse_EnsemblID)), c("Rat_EnsemblID", "Mouse_EnsemblID")]
convert_by_EnsemblID <- convert_by_EnsemblID[convert_by_EnsemblID$Rat_EnsemblID %in% Rat_EnsemblID, ]
mouse_name_map <- setNames(convert_by_EnsemblID$Mouse_EnsemblID, convert_by_EnsemblID$Rat_EnsemblID)

# Step 3: Map mouse names to all genes
mouse_names_for_all <- mouse_name_map[Rat_EnsemblID]  # NA if no mapping

# Step 4: Add as feature-level metadata
# Use Seurat::AddMetaData() on the assay
MS[["RNA"]] <- AddMetaData(
  object = MS[["RNA"]],
  metadata = data.frame(Mouse_EnsemblID = mouse_names_for_all, row.names = rownames(MS[["RNA"]])
  ))


# now make a new metadata column with mouse symbols
mouse_symbol_map <- convert[
  !is.na(convert$Mouse_EnsemblID) & !is.na(convert$Mouse_Symbol),
  c("Mouse_EnsemblID", "Mouse_Symbol")
]

mouse_symbol_map <- setNames(
  mouse_symbol_map$Mouse_Symbol,
  mouse_symbol_map$Mouse_EnsemblID
)

mouse_ens_ids <- MS[["RNA"]]@meta.data$Mouse_EnsemblID

mouse_symbols_for_all <- mouse_symbol_map[mouse_ens_ids]

MS[["RNA"]] <- AddMetaData(
  object = MS[["RNA"]],
  metadata = data.frame(
    Mouse_Symbol = mouse_symbols_for_all,
    row.names = rownames(MS[["RNA"]])
  )
)

# save MS8 with mouse genes in metadata -----------------------------------------
# saveRDS(MS, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS8.RDS")
# MS <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS8.RDS")
# -----------------------------------------------------------------------------


# add mapmycells class and subclass annotations to MS
unique(colnames(MS))

# Make a copy of MS column names for mapping
ms_cells_fixed <- colnames(MS)

# 1. Replace '-' with '.' for consistent formatting
 ms_cells_fixed <- gsub("-", ".", ms_cells_fixed)

# 2. Convert trailing ".1" if needed (adjust depending on your data)
 ms_cells_fixed <- gsub("_1$", ".1", ms_cells_fixed)
unique(rownames(mapmycells))

# Check intersection
cells_to_map <- intersect(ms_cells_fixed, rownames(mapmycells))
length(cells_to_map)  # should be >0 now

# Initialize metadata columns
MS$MapMyCells_Class <- NA
MS$MapMyCells_Subclass <- NA
MS$MapMyCells_Subclass_Probability <- NA

# Map annotations for matched cells
MS$MapMyCells_Class[match(cells_to_map, ms_cells_fixed)] <-
  mapmycells[cells_to_map, "class_name"]

MS$MapMyCells_Subclass[match(cells_to_map, ms_cells_fixed)] <-
  mapmycells[cells_to_map, "subclass_name"]

MS$MapMyCells_Subclass_Probability[match(cells_to_map, ms_cells_fixed)] <- 
  mapmycells[cells_to_map, "subclass_bootstrapping_probability"]


# save MS with mouse genes in metatdata and MapMyCells annotations -------------
# saveRDS(MS, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS9.RDS")
#  MS <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS9.RDS")
# ------------------------------------------------------------------------------

################################################################################
# Remove cells with uncertain cell types
################################################################################

# eliminate cells with MapMyCells_Subclass_Probability less than 95%

MS_filtered <- subset(
  MS,
  subset = MapMyCells_Subclass_Probability >= 0.95
)

num_cells <- nrow(MS_filtered@meta.data)
print(num_cells)

MS <- MS_filtered


# Count cells per cell type
cell_counts <- table(MS$MapMyCells_Subclass)

# Keep cell types with at least 10 cells
keep_types <- names(cell_counts[cell_counts >= 10])

# Subset Seurat object
MS <- subset(MS, subset = MapMyCells_Subclass %in% keep_types)

# 71086 cells after filtering



# save as MS11 -----------------------------------------------------------------
# saveRDS(MS, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS11.RDS")
# MS <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS11.RDS")
# ------------------------------------------------------------------------------

# make table of MapMyCells counts

library(gridExtra)
library(grid)


# Create table sorted descending
table <- as.data.frame(sort(table(MS$MapMyCells_Subclass), decreasing = TRUE))
colnames(table) <- c("Subclass", "n_cells")

# Create table grob with row heights proportional
table_grob <- tableGrob(
  table,
  rows = NULL,
  theme = ttheme_default(
    core = list(fg_params = list(cex = 0.8), # smaller font
                padding = unit(c(2,2), "mm")),
    colhead = list(fg_params = list(cex = 1))
  )
)

# Save directly to PDF, height scaled to number of rows
pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/celltypes.pdf", width = 6, height = nrow(table)*0.25 + 1)
grid.draw(table_grob)
dev.off()


# write_csv(table, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS95.csv")

# use csv to regroup cell types and eliminate/rename those that shouldn't be there (e.g. ctx)


# ------------------------------------------------------------------------------
# merge subtypes as defined in csv
# ------------------------------------------------------------------------------

library(dplyr)

subclass_to_myclass <- c(
  
  # Astrocytes
  "318 Astro-NT NN" = "Astrocytes",
  "319 Astro-TE NN" = "Astrocytes",
  
  # Cholinergic neurons
  "058 PAL-STR Gaba-Chol" = "Cholinergic Neurons",
  
  # GABAergic neurons
  "057 NDB-SI-MA-STRv Lhx8 Gaba" = "GABAergic Neurons (1)",
  "066 NDB-SI-ant Prdm12 Gaba"   = "GABAergic Neurons (2)",
  "067 LSX Sall3 Pax6 Gaba" = "GABAergic Neurons (2)",
  "060 OT D3 Folh1 Gaba" = "GABAergic Neurons (2)",
  "089 PVR Six3 Sox3 Gaba" = "GABAergic Neurons (2)",
  "064 STR-PAL Chst9 Gaba" = "GABAergic Neurons (2)",
  "090 BST-MPN Six3 Nrgn Gaba" = "GABAergic Neurons (2)",
  "069 LSX Nkx2-1 Gaba" = "GABAergic Neurons (2)",
  
  # Glutamatergic neurons
  "115 MS-SF Bsx Glut" = "Glutamatergic Neurons",
  "001 CLA-EPd-CTX Car3 Glut" = "Glutamatergic Neurons",
  "243 PGRN-PARN-MDRN Hoxb5 Glut" = "Glutamatergic Neurons",
  "119 SI-MA-LPO-LHA Skor1 Glut" = "Glutamatergic Neurons",
  
  # Immune
  "334 Microglia NN" = "Microglia",
  
  # Oligodendrocyte lineage
  "327 Oligo NN" = "Oligodendrocytes",
  "326 OPC NN"   = "OPCs",
  
  # Vascular
  "333 Endo NN" = "Endothelial Cells",
  "331 Peri NN" = "Vascular",
  "330 VLMC NN" = "Vascular"
)



MS$MyClass <- as.character(
  subclass_to_myclass[ MS$MapMyCells_Subclass ]
)


# Replace NA (unmapped subclasses) with "Others"
MS$MyClass[is.na(MS$MyClass)] <- "Others"

# remove others and remove NAs
MS <- subset(
  MS,
  subset = MyClass != "Others"
)


table(MS$MyClass, useNA = "ifany")

MS <- subset(
  MS,
  subset = !is.na(MyClass)
)

MS

# 70521 cells after filtering


################################################################################
# make new Immature Oligodendrocyte group in MyClass based on Bcas1 expression
################################################################################

Idents(MS) <- MS$MyClass
oligos <- subset(MS, idents = "Oligodendrocytes")

Bcas1_expr <- FetchData(oligos, vars = "Bcas1")

cutoff <- quantile(Bcas1_expr$Bcas1, 0.5)
Bcas1_pos <- Bcas1_expr$Bcas1 >= cutoff


oligos$MyClass_Bcas1 <- ifelse(
  Bcas1_pos,
  "Immature Oligodendrocytes",
  "Mature Oligodendrocytes"
)

MS$MyClass_Bcas1 <- MS$MyClass

MS$MyClass_Bcas1[colnames(oligos)] <- oligos$MyClass_Bcas1

Idents(MS) <- MS$MyClass_Bcas1

unique(MS$MyClass_Bcas1)

library(gridExtra)
library(grid)


# Create table sorted descending
table <- as.data.frame(sort(table(MS$MyClass_Bcas1), decreasing = TRUE))
colnames(table) <- c("MyClass_Bcas1", "n_cells")

# Create table grob with row heights proportional
table_grob <- tableGrob(
  table,
  rows = NULL,
  theme = ttheme_default(
    core = list(fg_params = list(cex = 0.8), # smaller font
                padding = unit(c(2,2), "mm")),
    colhead = list(fg_params = list(cex = 1))
  )
)

# Save directly to PDF, height scaled to number of rows
pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/celltypes.pdf", width = 6, height = nrow(table)*0.25 + 1)
grid.draw(table_grob)
dev.off()

# Save as MS13 -----------------------------------------------------------------
# saveRDS(MS, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS13.RDS")
# MS <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MS13.RDS")
# ------------------------------------------------------------------------------



################################################################################
# Cluster filtered MS and perform cell type verification
################################################################################

MS <- FindClusters(MS, reduction = "harmony")

Idents(MS) <- "group"
p <- DimPlot(MS,reduction = "umap",  label = FALSE)
p
ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/cluster.pdf",
  plot = p,
  height = 10,
  width = 12,
  dpi = 300
)


# make PCA plot

sample_to_group <- MS@meta.data |>
  dplyr::select(sample, group) |>
  dplyr::distinct() |>
  dplyr::arrange(group, sample)

sample_order <- sample_to_group$sample

p <- DotPlot(MS,
             features = gene_order,
             group.by = "sample",
             cols = c("gray90", "red"),
             col.min = 0,
             col.max = 2,
             dot.scale = 24)

sample_labels <- setNames(sample_to_group$group, sample_to_group$sample)

p$data$pct.exp[p$data$pct.exp < 0.1] <- NA

p <- p +
  scale_y_discrete(limits = sample_order, labels = sample_labels) +
  scale_size_area(max_size = 24, breaks = c(1,10,20,30))

p


# make dotplot

Idents(MS) <- MS$MyClass_Bcas1

DimPlot(MS, reduction = "pca", label = TRUE)

gene_order <-  c('Gfap', # 318 Astro-NT NN
                 
                 'Ngfr',  # 058 PAL-STR Gaba-Chol
                 
                 'Gad1', # Gaba
                 
                 'Slc17a7', 'Slc17a6', # 115 MS-SF Bsx Glut
                 
                 'P2ry12', # 334 Microglia NN
                 
                 'Cldn11', # 327 Oligo NN
                 
                 'Gpr17', 'Bcas1', # Immature Oligo
                 
                 'Cspg4', # 326 OPC NN
                 
                 'Cd34', # 333 Endo NN
                 
                 'Kcnj8', 'Ogn')# Vascular
                 
  
  

p <- DotPlot(MS, features = gene_order, group.by = "MyClass_Bcas1", cols = c("gray90", "red"),col.min = 0,col.max = 2,dot.scale=30) +
  scale_x_discrete(limits = gene_order) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 50),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.text  = element_text(size = 50),
        legend.title = element_text(size = 50, margin = margin(b = 30)),
        legend.spacing.y = unit(3, "cm") ) +
  guides(color = guide_colorbar(barheight = unit(10, "cm")))


p
ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/markersMS.pdf",
  plot = p,
  width = 30,
  height = 24,
  dpi = 300
)


# ------------------------------------------------------------------------------


################################################################################
# Construct Metacells via metacell
################################################################################

# BiocManager::install("tanaylab/metacell")
# BiocManager::install(c("impute", "preprocessCore", "GO.db"))
# install.packages(c('WGCNA', 'UCell', 'GeneOverlap'))
# install.packages("curl", repos = "https://jeroen.r-universe.dev")
# conda install conda-forge::r-enrichr
# BiocManager::install(c("WGCNA", "UCell", "GenomicRanges", "GeneOverlap"))
# devtools::install_github('smorabit/hdWGCNA', ref='dev',force = TRUE)

library(metacell)
library(Seurat)
library(hdWGCNA)
library(WGCNA)

# 1. Setting up the hdWGCNA experiment
# This step is the entry point for hdWGCNA; it adds an empty hdWGCNA experiment to your Seurat object.
seurat_obj <- SetupForWGCNA(
  MS,
  gene_select = "variable", # or "fraction", "topN", "variable"
  #  fraction = 0.05,          # If gene_select = "fraction", select the percentage 
  # of expressed genes. If your goal is to cluster and DEG only on major highly 
  # expressed genes, and the cell number is large, 5% can be used as an initial attempt.
  #To achieve more comprehensive coverage of subgroup markers, it is recommended 
  # to increase the fraction (e.g., 0.10~0.20), or directly use "variable" 
  # screening (e.g., for hypervariable genes).
  # GZ used "variable"
  wgcna_name = "Metacell" # Give your WGCNA experiment a name
)

# variable has 2000 genes.
length(seurat_obj@misc$Metacell$wgcna_genes)

parallel::detectCores() #check cpu core number

# optionally enable multithreading
enableWGCNAThreads(nThreads = parallel::detectCores())

# 2. Constructing Metacells
# Using MetacellsByGroups()
# If you want to construct Metacells within predefined cell groups (such as Seurat clusters, cell types),
# you can use MetacellsByGroups. This helps reduce sparsity while maintaining cell type specificity.
# For example, constructing metacells within each seurat_clusters:
seurat_obj <- MetacellsByGroups(
  seurat_obj,
  group.by = c("MyClass_Bcas1", "sample"), # Group by which metadata column 
  # Specify the columns in seurat_obj@meta.data to group by
  k = 10,                       # Number of nearest neighbor cells aggregated within each group's Metacell
  reduction = "pca",
  dims = 1:40,
  layer = "counts",
  mode = "sum",
  min_cells = 10,              # Minimum number of cells required to build a metacell
  max_shared = 1,
  target_metacells = 7052,    # set for ~ 10 cells per metacell
  ident.group ="MyClass_Bcas1", # set the Idents of the metacell seurat object
  verbose = TRUE
)



table(seurat_obj$MyClass_Bcas1, seurat_obj$group) # Observe the number of cells in each group

 
# The most convenient way is to directly extract the Seurat object of the 
# Metacell from the hdWGCNA experiment.

metacell_obj <- GetMetacellObject(seurat_obj)

# check the marker expression using metacell

Idents(metacell_obj) <- metacell_obj$MyClass_Bcas1

gene_order <-  c('Gfap', # 318 Astro-NT NN
                 
                 'Ngfr',  # 058 PAL-STR Gaba-Chol
                 
                 'Gad1', # Gaba
                 
                 'Slc17a7', 'Slc17a6', # 115 MS-SF Bsx Glut
                 
                 'P2ry12', # 334 Microglia NN
                 
                 'Cldn11', # 327 Oligo NN
                 
                 'Gpr17', 'Bcas1', # Immature Oligo
                 
                 'Cspg4', # 326 OPC NN
                 
                 'Cd34', # 333 Endo NN
                 
                 'Kcnj8', 'Ogn')# Vascular



p <- DotPlot(metacell_obj, features = gene_order,cols = c("gray90", "red"),col.min = 0,col.max = 2,dot.scale=10) +
  scale_x_discrete(limits = gene_order) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.text  = element_text(size = 24),
        legend.title = element_text(size = 24))
p
ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/markers.pdf",
  plot = p,
  width = 16,
  height = 10,
  dpi = 300
)

library(dplyr)
library(gridExtra)

# Metacell counts
mc_counts <- as.data.frame(table(metacell_obj$MyClass_Bcas1))
colnames(mc_counts) <- c("MyClass_Bcas1", "n_metacells")

# Single-cell counts (from Seurat object)
sc_counts <- MS@meta.data %>%
  dplyr::group_by(MyClass_Bcas1) %>%
  dplyr::summarise(n_cells = n(), .groups = "drop")

# Merge + sort
table <- mc_counts %>%
  left_join(sc_counts, by = "MyClass_Bcas1") %>%
  arrange(desc(n_metacells))

# Create table grob
table_grob <- tableGrob(
  table,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(cex = 0.8),
      padding = unit(c(2, 2), "mm")
    ),
    colhead = list(fg_params = list(cex = 1))
  )
)


# Save directly to PDF, height scaled to number of rows
pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/celltypes.pdf", width = 6, height = nrow(table)*0.25 + 1)
grid.draw(table_grob)
dev.off()

metacell_obj

ls("package:metacell", pattern="*Metacell_Object*")

### save seurat_object

library(Seurat)
library(hdWGCNA)
library(WGCNA)

# The most convenient way is to directly extract the Seurat object of the Metacell from the hdWGCNA experiment.
metacell_obj <- GetMetacellObject(seurat_obj)

library(stringr)
replace_vec <- c('Lane-1__MS1__115-137'='XYT', 'Lane-1__MS2__119-152'='XXT', 'Lane-1__MS3__173-235'='XXO',
              'Lane-1__MS4__149-159'='XYT', 'Lane-1__MS5__171-129'='XXT', 'Lane-1__MS6__191-174'='XXO',
              'Lane-1__MS7__186-197'='XYT', 'Lane-1__MS8__123-165'='XXT', 'Lane-1__MS9__190-194'='XXO')

metacell_obj$group <- str_replace_all(metacell_obj$sample, replace_vec)

library(stringr)
replace_vec <- c('Lane-1__MS1__115-137'='batch1', 'Lane-1__MS2__119-152'='batch1', 'Lane-1__MS3__173-235'='batch1',
                 'Lane-1__MS4__149-159'='batch2', 'Lane-1__MS5__171-129'='batch2', 'Lane-1__MS6__191-174'='batch2',
                 'Lane-1__MS7__186-197'='batch3', 'Lane-1__MS8__123-165'='batch3', 'Lane-1__MS9__190-194'='batch3')

metacell_obj$batch <- str_replace_all(metacell_obj$sample, replace_vec) 


# transfer mouse symbol metadata from MS to metacell_obj
mouse_symbol <- MS@assays[["RNA"]]@meta.data$Mouse_Symbol
names(mouse_symbol) <- rownames(MS[["RNA"]])

stopifnot(
  identical(
    rownames(MS[["RNA"]]),
    rownames(metacell_obj[["RNA"]])
  )
)

metacell_obj[["RNA"]] <- AddMetaData(
  object = metacell_obj[["RNA"]],
  metadata = data.frame(
    Mouse_Symbol = mouse_symbol,
    row.names = rownames(metacell_obj[["RNA"]])
  )
)


# use the TransferData() from seurat for cell annotation
MSMC <- NormalizeData(metacell_obj)
MSMC <- FindVariableFeatures(MSMC)
MSMC <- ScaleData(MSMC)
MSMC <- RunPCA(MSMC) #default 50 PC
MSMC <- RunHarmony(MSMC, group.by.vars='batch')
MSMC <- FindNeighbors(MSMC, reduction = "harmony", dims = 1:40)
MSMC <- FindClusters(MSMC)
MSMC <- RunUMAP(MSMC, dims = 1:40)


Idents(MSMC) <- "MyClass_Bcas1"

Idents(MSMC) <- factor(Idents(MSMC))
cluster_names <- levels(Idents(MSMC))
n <- length(cluster_names)

umap <- Embeddings(MSMC, "umap")

centroids <- aggregate(
  umap,
  by = list(cluster = Idents(MSMC)),
  FUN = mean
)

ordered_clusters <- centroids$cluster
n <- length(ordered_clusters)

base_cols <- scales::hue_pal()(n)

reorder_idx <- c(seq(1, n, by = 2), seq(2, n, by = 2))

cols <- setNames(
  base_cols[reorder_idx],
  ordered_clusters
)
p <- DimPlot(MSMC,reduction = "umap",  label = TRUE, cols = cols)

p

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/cluster.pdf",
  plot = p,
  width = 9,
  dpi = 300
)


Idents(MSMC) <- "batch"
p <- DimPlot(MSMC,reduction = "umap",  label = FALSE)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/batch.pdf",
  plot = p,
  dpi = 300
)

Idents(MSMC) <- "group"
p <- DimPlot(MSMC,reduction = "umap",  label = FALSE)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/group.pdf",
  plot = p,
  dpi = 300
)


# Save as MSMC -----------------------------------------------------------------
# saveRDS(MSMC, file = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MSMC.RDS")
# MSMC <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MSMC.RDS")
# ------------------------------------------------------------------------------


################################################################################
#  DEG Analysis with Limma Trend 
################################################################################

##Script to run limma trend and limma voom for DEG analysis

library(limma)
library(edgeR)
library(data.table)
library(dplyr)
library(plyr)
library(Seurat)
library(org.Mm.eg.db)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggpubr)
library(DESeq2)


##load seurat metacell object 


celltypes <- unique(MSMC$MyClass_Bcas1)

deg_list <- list()

assay_to_use <- "RNA"

for (ct in celltypes) {
  
  message("Processing MyClass_Bcas1: ", ct)
  
  ## Subset by cell type and genotype
  obj <- subset(
    MSMC,
    subset = MyClass_Bcas1 == ct & group %in% c("XXO", "XXT", "XYT")
  )
  
  ## Skip if too few cells
  if (ncol(obj) < 10) {
    message("  Skipping: too few cells")
    next
  }
  
  meta <- obj@meta.data
  meta$group <- factor(meta$group, levels = c("XXO", "XXT", "XYT"))
  meta$batch <- factor(meta$batch)
  
  ## Must have all 3 genotypes for pairwise contrasts
  if (nlevels(droplevels(meta$group)) < 2) {
    message("  Skipping: missing genotypes")
    next
  }
  
  ## Get counts (Seurat v5 compatible)
  gene_counts <- GetAssayData(
    obj,
    assay = assay_to_use,
    slot  = "counts"
  )
  gene_counts <- round(gene_counts)
  
  ## edgeR filtering: keep genes expressed (CPM >= 1) in at least 2 samples
  dge0 <- DGEList(gene_counts)
  cpm_mat <- edgeR::cpm(dge0)
  keep <- apply(cpm_mat, 1, max) >= 2
  dge <- dge0[keep, , keep.lib.sizes = FALSE]
  
  ## ---- Labeled output ----
  cat("\n===== DEG FILTERING SUMMARY =====\n")
  cat("Total genes in assay (before filtering): ",
      nrow(dge0), "\n")
  
  cat("Genes expressed at CPM ≥ 2 in ≥1 sample (tested genes): ",
      nrow(dge), "\n")
  
  cat("Genes removed (low / no expression): ",
      nrow(dge0) - nrow(dge), "\n")
  
  cat("Percent of genes retained for DEG testing: ",
      round(100 * nrow(dge) / nrow(dge0), 2), "%\n")
  
  cat("=================================\n\n")
  
  if (nrow(dge) < 10) {
    message("  Skipping: too few genes after filtering")
    next
  }
  
  dge <- calcNormFactors(dge)
  
  ## Design matrix (include batch only if usable)
  if (nlevels(meta$batch) > 1) {
    design <- model.matrix(~ 0 + group + batch, data = meta)
  } else {
    design <- model.matrix(~ 0 + group, data = meta)
  }
  
  ## Expression matrix for limma
  y <- new("EList")
  y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
  
  fit <- lmFit(y, design)
  
  ## Pairwise contrasts
  contrasts <- makeContrasts(
    XYT_vs_XXO = groupXYT - groupXXO,
    XYT_vs_XXT = groupXYT - groupXXT,
    XXT_vs_XXO = groupXXT - groupXXO,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  
  ## Collect results
  resultlist <- list(
    XYT_vs_XXO = topTable(fit2, coef = "XYT_vs_XXO", n = Inf, adjust.method = "BH"),
    XYT_vs_XXT = topTable(fit2, coef = "XYT_vs_XXT", n = Inf, adjust.method = "BH"),
    XXT_vs_XXO = topTable(fit2, coef = "XXT_vs_XXO", n = Inf, adjust.method = "BH")
  )
  
  deg_list[[ct]] <- resultlist
  
  ## Save per cell type
  save(
    resultlist,
    file = paste0(
      "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/",
      gsub(" ", "_", ct),
      "_pairwise_genotype.rda"
    )
  )
}



saveRDS(deg_list, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/Downsample_Limma_batch.RDS")



# ------------------------------------------------------------------------------

library(dplyr)

# Assume the deg_list structure is deg_list[[ct]][[comp]], and each df row name 
# is either the gene name or the Ensembl ID.
deg_all <- do.call(
  rbind,
  lapply(names(deg_list), function(ct) {
    lapply(names(deg_list[[ct]]), function(comp) {
      df <- deg_list[[ct]][[comp]]
      df$gene <- rownames(df)      # Put the gene name in a separate column
      df$CellType <- ct
      df$Comparison <- comp
      df
    })
  }) %>% unlist(recursive = FALSE)
)

# Row names are currently numbers, and gene names are in the gene column. 
# You can use dplyr::select to adjust the order.
deg_all <- deg_all %>% dplyr::select(gene, CellType, Comparison, everything())

write.csv(deg_all, file = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/all_DEG_results_Limma.csv", row.names = FALSE)

# plot DEG number

library(readxl)
library(dplyr)
library(ggplot2)

# 2. Set the file path

# Please replace "your_deg_file.xlsx" with your actual Excel file name and path.
file_path <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/all_DEG_results_Limma.csv"

# 3. Read Excel file
deg_data <- read.csv(file_path)

# 4. Filter for significant DEGs and count the number of up-regulations/down-regulations

# Define significance threshold and logFC threshold
p_val_adj_threshold <- 0.05    # Adjusted p-value threshold
log2FC_threshold <- 0.1       # avg_log2FC Absolute value greater than 0.1

# Perform filtering and statistics
deg_summary_filtered <- deg_data %>%
  dplyr::filter(adj.P.Val < p_val_adj_threshold) %>%
  dplyr::mutate(
    regulation = dplyr::case_when(
      logFC >  log2FC_threshold  ~ "Up-regulated",
      logFC < -log2FC_threshold  ~ "Down-regulated",
      TRUE                       ~ "Not significant change"
    )
  ) %>%
  dplyr::filter(regulation != "Not significant change") %>%
  dplyr::group_by(CellType, regulation) %>%
  dplyr::summarise(
    count = dplyr::n(),
    .groups = "drop"
  )


# 5. View Results
print(paste0("Statistical results (p_val_adj < ", p_val_adj_threshold, " and |avg_log2FC| > ", log2FC_threshold, "):"))
print(deg_summary_filtered)

# --- Define Color Mapping ---

# Ensure the names here are exactly the same as those you created in the regulation column.

# Here, we'll set Up-regulated to blue and Down-regulated to red, as per your request.

color_mapping <- c("Up-regulated" = "red", "Down-regulated" = "blue")

# If you prefer the standard: Up-regulated red, Down-regulated blue, you can set it like this:

# color_mapping <- c("Up-regulated" = "red", "Down-regulated" = "blue")
# --------------------

# 6. Visualize the results (optional)

# Create a stacked bar chart
p <- ggplot(deg_summary_filtered, aes(x = CellType, y = count, fill = regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = paste0("Number of Up-regulated and Down-regulated DEGs per Cell Type\n(p_val_adj < ", p_val_adj_threshold, " & |avg_log2FC| > ", log2FC_threshold, ")"),
    x = "Cell Type",
    y = "Number of DEGs",
    fill = "Regulation"
  ) +
  scale_y_continuous(expand = c(0, 0)) + # Ensure the y-axis starts from 0 and there is no extra blank space at the top.
  scale_x_discrete(expand = c(0, 0)) +   # Eliminate blank spaces on both sides of the x-axis
  scale_fill_manual(values = color_mapping) + # <-- **Key Modification: Manually Set Colors**
  theme_minimal() + # Basic Theme
  theme(
    panel.grid.major = element_blank(), # Remove main grid lines
    panel.grid.minor = element_blank(), # Remove secondary grid lines
    axis.line = element_line(colour = "black"), # Add coordinate axes and set them to black
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), # Rotate the x-axis labels and set their color.
    axis.text.y = element_text(color = "black"), # Set y-axis label color
    axis.title.x = element_text(color = "black"), # Set x-axis title color
    axis.title.y = element_text(color = "black")  # Set y-axis title color
  )

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/Limma_DEGs.pdf",
  plot = p,
  width = 6, height = 4, units = "in",
  dpi = 300
)



################################################################################
# # plot DEG related figures
################################################################################

# Set file path
file_path <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/all_DEG_results_Limma.csv"

# 3. Read Excel file
deg_data <- read.csv(file_path)

deg_list <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/Downsample_Limma_batch.RDS")

library(dplyr)

# Assume the deg_list structure is deg_list[[ct]][[comp]], and each df row name 
# is either the gene name or the Ensembl ID.
deg_all <- do.call(
  rbind,
  lapply(names(deg_list), function(ct) {
    lapply(names(deg_list[[ct]]), function(comp) {
      df <- deg_list[[ct]][[comp]]
      df$gene <- rownames(df)     # Put the gene name in a separate column
      df$CellType <- ct
      df$Comparison <- comp
      df
    })
  }) %>% unlist(recursive = FALSE)
)

# Row names are currently numbers, and gene names are in the gene column. 
# You can use dplyr::select to adjust the order.
deg_all <- deg_all %>% dplyr::select(gene, CellType, Comparison, everything())

rownames(deg_all) <- NULL

# two types of gene names
deg_all <- deg_all %>%
  mutate(gene_type = ifelse(grepl("^ENSRNOG", gene), "Ensembl", "Symbol"))

ens_ids <- unique(deg_all$gene[deg_all$gene_type == "Ensembl"])
symbols <- unique(deg_all$gene[deg_all$gene_type == "Symbol"])

library(biomaRt)
rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")

# Use Ensembl ID to search
ens_annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
  filters = "ensembl_gene_id",
  values = ens_ids,
  mart = rat
)

# Use symbol to search
symbol_annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
  filters = "external_gene_name",
  values = symbols,
  mart = rat
)

# combine two tables
gene_annot <- dplyr::bind_rows(
  ens_annot %>% dplyr::rename(gene = ensembl_gene_id),
  symbol_annot %>% dplyr::rename(gene = external_gene_name)
)

deg_all2 <- deg_all %>%
  left_join(gene_annot, by = "gene") %>%
  mutate(class = ifelse(chromosome_name %in% c("X", "Y"), "Sex-linked", "Autosomal"))

# save deg_all2 ----------------------------------------------------------------
# saveRDS(deg_all2, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg_all2.RDS")
# deg_all2 <- readRDS("/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Limma/deg_all2.RDS")
# ------------------------------------------------------------------------------

write.csv(deg_all2, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEGs.csv")

library(dplyr)
library(ggplot2)
library(forcats)

# 1. Keep only the group XYT_vs_XXO
deg_plot <- deg_all2 %>%
  filter(Comparison == "XYT_vs_XXT" & adj.P.Val < 0.05)

# 2. Group statistics by CellType and class
plotdat <- deg_plot %>%
  dplyr::group_by(CellType, class) %>%
  dplyr::summarise(Freq = n()) %>%
  ungroup()

# 3. Sort CellType in descending order of Freq.
plotdat <- plotdat %>%
  mutate(CellType = fct_reorder(CellType, Freq, .desc = TRUE))
plotdat$facet_label <- "SCC (XYT vs XXT)"

# 4. Draw the drawing and adjust the font size.
p <- ggplot(plotdat, aes(x = CellType, y = log2(Freq), fill = class)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "log2(DEG Count)", fill = "Cell Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), # Adjust font size
        axis.text.y = element_text(angle = 0, hjust = 1, size = 10),
        panel.grid.major = element_blank(),   # Remove main grid lines
        panel.grid.minor = element_blank()  ) + # Remove secondary grid lines
  facet_wrap(~facet_label)

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_plot.pdf",
  plot = p,
  width = 6, height = 4, units = "in",
  dpi = 300
)

#-------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(forcats)

# 1. Keep only the XYT_vs_XXO group, padj < 0.05
deg_plot <- deg_all2 %>%
  filter(Comparison == "XYT_vs_XXO" & adj.P.Val < 0.05)

# 2. Add topping/down tags
deg_plot <- deg_plot %>%
  mutate(Direction = ifelse(logFC > 0, "Up", "Down"))

# 3. Count the number of genes upregulated and downregulated for each cell type
deg_plot <- deg_plot %>%
  dplyr::mutate(Direction = ifelse(logFC > 0, "Up", "Down"))

plotdat <- deg_plot %>%
  dplyr::group_by(CellType, Direction) %>%
  dplyr::summarise(Freq = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(Freq = ifelse(Freq == 0, NA, Freq)) %>%   # avoid log2(0)
  tidyr::pivot_wider(
    names_from = Direction,
    values_from = Freq,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    Up_log2   = ifelse(Up   > 0,  log2(Up),        NA),
    Down_log2 = ifelse(Down > 0, -log2(Down),      NA),
    Total     = Up + Down
  ) %>%
  tidyr::pivot_longer(
    cols = c("Up_log2", "Down_log2"),
    names_to = "Direction",
    values_to = "log2Freq"
  ) %>%
  dplyr::mutate(
    Direction = dplyr::recode(Direction,
                              "Up_log2" = "Up",
                              "Down_log2" = "Down"),
    CellType  = forcats::fct_reorder(CellType, abs(Total), .desc = TRUE),
    facet_label = "SCC (XYT vs XXT)"
  ) %>%
  dplyr::filter(!is.na(log2Freq))

#2. Drawing
p <- ggplot(plotdat, aes(x = CellType, y = log2Freq, fill = Direction)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "log2(DEG Count)", x = NULL, fill = "Direction") +
  facet_wrap(~facet_label) +
  scale_fill_manual(values = c("Down" = "#00bfc4", "Up" = "#f8766d")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_plot6.pdf",
  plot = p,
  width = 6, height = 4, units = "in",
  dpi = 300
)

# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)

# 1. Keep only the XXT_vs_XXO group
deg_sub <- deg_all2 %>% filter(Comparison == "XXT_vs_XXO")

# 2. Extract autosomal and sex-linked data separately.
deg_autosomal <- deg_sub %>% filter(class == "Autosomal")
deg_sexlinked <- deg_sub %>% filter(class == "Sex-linked")

# 3. Aggregate logFC for each (gene, CellType) (take the mean to prevent duplication)
deg_autosomal_unique <- deg_autosomal %>%
  dplyr::group_by(gene, CellType) %>%
  dplyr::summarise(
    logFC = mean(logFC, na.rm = TRUE),
    .groups = "drop"
  )

deg_sexlinked_unique <- deg_sexlinked %>%
  dplyr::group_by(gene, CellType) %>%
  dplyr::summarise(logFC = mean(logFC, na.rm = TRUE), .groups = 'drop')

# 4. Wide table transformation
logFC_mat_autosomal <- deg_autosomal_unique %>%
  pivot_wider(names_from = CellType, values_from = logFC)

logFC_mat_sexlinked <- deg_sexlinked_unique %>%
  pivot_wider(names_from = CellType, values_from = logFC)

# 5. Set row names to numeric matrix
logFC_mat_autosomal <- as.data.frame(logFC_mat_autosomal)
rownames(logFC_mat_autosomal) <- logFC_mat_autosomal$gene
logFC_mat_autosomal$gene <- NULL
logFC_mat_autosomal[] <- lapply(logFC_mat_autosomal, as.numeric)
logFC_mat_autosomal <- as.matrix(logFC_mat_autosomal)

logFC_mat_sexlinked <- as.data.frame(logFC_mat_sexlinked)
rownames(logFC_mat_sexlinked) <- logFC_mat_sexlinked$gene
logFC_mat_sexlinked$gene <- NULL
logFC_mat_sexlinked[] <- lapply(logFC_mat_sexlinked, as.numeric)
logFC_mat_sexlinked <- as.matrix(logFC_mat_sexlinked)

# 6. Correlation Calculation
cor_mat_autosomal <- cor(logFC_mat_autosomal, use = "pairwise.complete.obs", method = "pearson")
cor_mat_sexlinked <- cor(logFC_mat_sexlinked, use = "pairwise.complete.obs", method = "pearson")

cor_mat_autosomal

get_cor_pmat <- function(mat) {
  n <- ncol(mat)
  pmat <- matrix(NA, n, n)
  colnames(pmat) <- rownames(pmat) <- colnames(mat)
  for(i in 1:n) {
    for(j in 1:n) {
      if(i == j) {
        pmat[i, j] <- 1
      } else {
        test <- cor.test(mat[,i], mat[,j], method = "pearson")
        pmat[i, j] <- test$p.value
      }
    }
  }
  pmat
}
p_mat_autosomal <- get_cor_pmat(logFC_mat_autosomal)
p_mat_sexlinked <- get_cor_pmat(logFC_mat_sexlinked)

p_mat_autosomal

library(pheatmap)

# Prepare the asterisk annotation matrix
get_star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}
ann_autosomal <- get_star(p_mat_autosomal)
ann_sexlinked <- get_star(p_mat_sexlinked)

annotation_col = data.frame(
  CellType = colnames(cor_mat_autosomal),
  Effect = effect_vec[colnames(cor_mat_autosomal)]
)
annotation_row = data.frame(
  CellType = rownames(cor_mat_autosomal),
  Effect = effect_vec[rownames(cor_mat_autosomal)]
)

rownames(annotation_col) <- colnames(cor_mat_autosomal)
rownames(annotation_row) <- rownames(cor_mat_autosomal)

unique(annotation_col$CellType)


ann_colors <- list(
  CellType =  c( "Mature Oligodendrocytes" = "darkred",
                 "Immature Oligodendrocytes" = "red",
                 "GABAergic Neurons (1)" = "purple",
                 "GABAergic Neurons (2)" = "lavender",
                 "Astrocytes" = "yellow",
                 "OPCs" = "orange",
                 "Glutamatergic Neurons" = "green",
                 "Cholinergic Neurons" = "blue"
  ), # Your CellType color chart
  Effect = c("Normative" = "purple", "Other" = "lavender") # Custom
)

# heatmap will fail. ignore errors and move on

p <- pheatmap(
  cor_mat_autosomal,
  display_numbers = ann_autosomal,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  main = "Typical (XYT vs XXO) logFC Correlations (Autosomal)",
  fontsize_number = 12,
  breaks = seq(-1, 1, length.out = 101)
)


library(ComplexHeatmap)


pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/Correlation_autosomal.pdf", width=14, height=6) ## ggsave can't save not ggplot
print(p)
dev.off()

# -----------------------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)

# Assuming ann_colors, annotation_col, annotation_row, cor_mat_autosomal,
# and ann_autosomal are already prepared.

ht <- Heatmap(
  cor_mat_autosomal,
  name = "Correlation",  # Main legend title
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  top_annotation = HeatmapAnnotation(
    CellType = annotation_col$CellType,
    col = ann_colors,
    annotation_legend_param = list(
      CellType = list(title = "CellType")
    )
  ),
  left_annotation = rowAnnotation(
    CellType = annotation_row$CellType,
    Effect = annotation_row$Effect,
    col = ann_colors,
    annotation_legend_param = list(
      CellType = list(title = "CellType")
    )
  ),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Gonadal (XXT vs XXO) logFC Correlations
  (Autosomal)",
  heatmap_legend_param = list(
    title = "Correlation",
    at = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    legend_height = unit(2, "cm"),
    border = "black"
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(ann_autosomal[i, j], x, y, gp = gpar(fontsize = 10))
  }
)

# All legends are listed in the column on the right.
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE)

pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/Correlation_autosomal.pdf", width=8, height=7) ## ggsave can't save not ggplot
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

# ------------------------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)
library(grid) # for unit()

# HeatmapAnnotation and rowAnnotation only need to be defined once (consistent with autosomal).）
ht_sexlinked <- Heatmap(
  cor_mat_sexlinked,
  name = "Correlation",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  top_annotation = HeatmapAnnotation(
    CellType = annotation_col$CellType,
    col = ann_colors,
    annotation_legend_param = list(
      CellType = list(title = "CellType")
    )
  ),
  left_annotation = rowAnnotation(
    CellType = annotation_row$CellType,
    Effect = annotation_row$Effect,
    col = ann_colors,
    annotation_legend_param = list(
      CellType = list(title = "CellType")
    )
  ),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Gonadal (XXT vs XXO) logFC Correlations
  (Sex-linked)",
  heatmap_legend_param = list(
    title = "Correlation",
    at = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    legend_height = unit(2, "cm"), # Legend length can be adjusted as needed
    border = "black"
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(ann_sexlinked[i, j], x, y, gp = gpar(fontsize = 10))
  }
)

draw(ht_sexlinked, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE)

pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/Correlation_sex.pdf", width=8, height=7) ## ggsave can't save not ggplot
draw(ht_sexlinked, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()


# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(forcats)

# 1. Only keep 3 comparisons & FDR < 0.05
deg_plot <- deg_all2 %>%
  filter(
    Comparison %in% c("XYT_vs_XXO", "XXT_vs_XXO", "XYT_vs_XXT"),
    adj.P.Val < 0.05
  )

# 2. 设定每个comparison的facet标签，两行显示
comparison_labels <- c(
  "XYT_vs_XXO" = "Typical\n(XYT vs XXO)",
  "XXT_vs_XXO" = "Gonadal\n(XXT vs XXO)",
  "XYT_vs_XXT" = "SCC\n(XYT vs XXT)"
)
deg_plot$facet_label <- comparison_labels[deg_plot$Comparison]

# 3. Group statistics by CellType
plotdat <- deg_plot %>%
  dplyr::group_by(facet_label, CellType) %>%
  dplyr::summarise(Freq = n(), .groups = "drop")

# 4. Set the facet order (Normative, Gonad, Sex chromosome)
plotdat$facet_label <- factor(plotdat$facet_label, levels = c(
  "Typical\n(XYT vs XXO)",
  "Gonadal\n(XXT vs XXO)",
  "SCC\n(XYT vs XXT)"
))

# 5. Set CellType order in descending order of the number of DEG values for Normative features.
normative_order <- plotdat %>%
  filter(facet_label == "Typical\n(XYT vs XXO)") %>%
  arrange(desc(Freq)) %>%
  pull(CellType)
plotdat$CellType <- factor(plotdat$CellType, levels = normative_order)

# 6. Custom Colors
cell_colors <- c(
  "Mature Oligodendrocytes" = "darkred",
  "Immature Oligodendrocytes" = "red",
  "GABAergic Neurons (1)" = "purple",
  "GABAergic Neurons (2)" = "lavender",
  "Astrocytes" = "yellow",
  "OPCs" = "orange",
  "Glutamatergic Neurons" = "green",
  "Cholinergic Neurons" = "blue",
  "Endothelian Cells" = "lightblue",
  "Microglia" = "darkgreen",
  "Vascular" = "gray"
)

# 7. Drawing
p <- ggplot(plotdat, aes(x = CellType, y = log2(Freq), fill = CellType)) +
  geom_bar(stat = "identity", color = "black") +
  theme_bw() +
  labs(y = "log2(DEG Count)", x = "Cell Type", fill = "Cell Type") +
  scale_fill_manual(values = cell_colors, drop = FALSE) +
  scale_y_continuous(expand = c(0, 0), labels = as.integer) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  facet_wrap(~facet_label, nrow = 1, scales = "fixed")

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_plot8.pdf",
  plot = p,
  width = 8, height = 4, units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(forcats)

# 1. 
deg_plot <- deg_all2 %>%
  filter(
    Comparison %in% c("XYT_vs_XXO", "XXT_vs_XXO", "XYT_vs_XXT"),
    adj.P.Val < 0.05
  ) %>%
  mutate(Direction = ifelse(logFC > 0, "Up", "Down")) 

# 2. 
comparison_labels <- c(
  "XYT_vs_XXO" = "Typical\n(XYT vs XXO)",
  "XXT_vs_XXO" = "Gonadal\n(XXT vs XXO)",
  "XYT_vs_XXT" = "SCC\n(XYT vs XXT)"
)
deg_plot$facet_label <- comparison_labels[deg_plot$Comparison]

# 3. 
plotdat <- deg_plot %>%
  dplyr::group_by(facet_label, CellType, Direction) %>%
  dplyr::summarise(Freq = n(), .groups = "drop")

# 4. 
plotdat$facet_label <- factor(plotdat$facet_label, levels = c(
  "Typical\n(XYT vs XXO)",
  "Gonadal\n(XXT vs XXO)",
  "SCC\n(XYT vs XXT)"
))



# 5. 
normative_order <- plotdat %>%
  dplyr::filter(facet_label == "Typical\n(XYT vs XXO)") %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(
    Total = sum(Freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(Total)) %>%
  dplyr::pull(CellType)

plotdat$CellType <- factor(plotdat$CellType, levels = normative_order)

# 6.
plotdat <- plotdat %>%
  mutate(log2Freq = ifelse(Direction == "Down", -log2(Freq), log2(Freq)))

# 7. 
direction_colors <- c( "Up" = "#f8766d","Down" = "#00bfc4")  # Down=青, Up=粉


# 8. 
ggplot(plotdat, aes(x = CellType, y = log2Freq, fill = Direction)) +
  geom_bar(stat = "identity", color = "black") +
  theme_bw() +
  labs(y = "log2(DEG Count)", x = "Cell Type", fill = "Direction") +
  scale_fill_manual(values = direction_colors) +
  scale_y_continuous(expand = c(0, 0), labels = as.integer) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  facet_wrap(~facet_label, nrow = 1, scales = "fixed")

ggsave("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_all_group_direction.pdf",height = 6,width = 7)

# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(forcats)

# 1. 
deg_plot <- deg_all2 %>%
  filter(
    Comparison %in% c("XYT_vs_XXO", "XXT_vs_XXO", "XYT_vs_XXT"),
    adj.P.Val < 0.05
  )

# 2. 
comparison_labels <- c(
  "XYT_vs_XXO" = "Typical\n(XYT vs XXO)",
  "XXT_vs_XXO" = "Gonadal\n(XXT vs XXO)",
  "XYT_vs_XXT" = "SCC\n(XYT vs XXT)"
)
deg_plot$facet_label <- comparison_labels[deg_plot$Comparison]

# 3. 
plotdat <- deg_plot %>%
  dplyr::group_by(facet_label, CellType, class) %>%
  dplyr::summarise(Freq = n(), .groups = "drop")

# 4.
plotdat$facet_label <- factor(plotdat$facet_label, levels = c(
  "Typical\n(XYT vs XXO)",
  "Gonadal\n(XXT vs XXO)",
  "SCC\n(XYT vs XXT)"
))

# 5. 
normative_order <- plotdat %>%
  dplyr::filter(facet_label == "Typical\n(XYT vs XXO)") %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(Total = sum(Freq)) %>%
  dplyr::arrange(desc(Total)) %>%
  dplyr::pull(CellType)
plotdat$CellType <- factor(plotdat$CellType, levels = normative_order)

# 6. 
plotdat <- plotdat %>%
  mutate(log2Freq = log2(Freq))

# 7.
cell_class_colors <- c("Autosomal" = "#FF6F6F", "Sex-linked" = "#10CFC9")

# 8. 
ggplot(plotdat, aes(x = CellType, y = log2Freq, fill = class)) +
  geom_bar(stat = "identity", color = "black", position = "stack") +
  theme_bw() +
  labs(y = "log2(DEG Count)", x = "Cell Type", fill = "Cell Class") +
  scale_fill_manual(values = cell_class_colors) +
  scale_y_continuous(expand = c(0, 0), labels = as.integer) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  facet_wrap(~facet_label, nrow = 1, scales = "fixed")

ggsave("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_all_group_sexauto.pdf",height = 6,width = 7)

# ------------------------------------------------------------------------------

# # DEG across factors

library(dplyr)
deg_sig <- deg_all2 %>% 
  filter(adj.P.Val < 0.05) %>% 
  dplyr::select(gene, CellType, Comparison,class) # need to add dplyr::, otherwise error

library(dplyr)

dup_counts <- deg_sig %>%
  group_by(gene, CellType, Comparison,class) %>%
  tally() %>%
  filter(n > 1)

# View duplicated combinations and how many times they appear
print(dup_counts)
# cause of your earlier pivot_wider error.tidyr::pivot_wider()`: ! 
# Can't convert `fill` <logical> to <list>. Duplicates must be resolved before
# reshaping the data.

library(dplyr)

deg_sig_nodup <- deg_sig %>%
  distinct(gene, CellType, Comparison,class, .keep_all = TRUE)

sum(duplicated(deg_sig_nodup[, c("gene", "CellType", "Comparison","class")]))
# Should be 0

library(tidyr)
deg_wide <- deg_sig_nodup %>%
  mutate(flag = TRUE) %>%
  pivot_wider(names_from = Comparison, values_from = flag, values_fill = FALSE)

# define genes category
deg_wide <- deg_wide %>%
  mutate(
    deg_class = case_when(
      XYT_vs_XXO & !XXT_vs_XXO & !XYT_vs_XXT ~ "Typical_unique",
      !XYT_vs_XXO & XXT_vs_XXO & !XYT_vs_XXT ~ "Gonadal",
      !XYT_vs_XXO & !XXT_vs_XXO & XYT_vs_XXT ~ "SCC",
      (XYT_vs_XXO + XXT_vs_XXO + XYT_vs_XXT) > 1 ~ "Shared"
    )
  )

library(dplyr)


deg_stats <- deg_wide %>%
  dplyr::group_by(CellType, class, deg_class) %>%
  dplyr::summarise(gene_count = n(), .groups = "drop") %>%
  dplyr::group_by(CellType, class) %>%
  dplyr::mutate(prop = gene_count / sum(gene_count))

deg_stats <- deg_stats %>%
  mutate(deg_class = factor(deg_class, levels = c("Typical_unique", "Shared", "Gonadal", "SCC")))

deg_stats_auto <- deg_stats %>% filter(class == "Autosomal")
deg_stats_sex <- deg_stats %>% filter(class == "Sex-linked")

# Autosomal
order_auto <- deg_stats_auto %>%
  filter(deg_class == "Typical_unique") %>%
  arrange(prop) %>%
  pull(CellType)
# to avoid missing cell type in factors, add cell type without Normative_unique
all_celltypes <- unique(deg_stats_auto$CellType)
order_auto_full <- c(setdiff(all_celltypes, order_auto),order_auto)
deg_stats_auto <- deg_stats_auto %>%
  mutate(CellType = factor(CellType, levels = order_auto_full))

# Sex-linked
order_sex <- deg_stats_sex %>%
  filter(deg_class == "Typical_unique") %>%
  arrange(prop) %>%
  pull(CellType)
# to avoid missing cell type in factors, add cell type without Normative_unique

all_celltypes <- unique(deg_stats_sex$CellType)
order_sex_full <- c(setdiff(all_celltypes, order_sex),order_sex)
deg_stats_sex <- deg_stats_sex %>%
  mutate(CellType = factor(CellType, levels = order_sex_full))

order_auto <- sort(unique(as.character(deg_stats_auto$CellType)))
order_sex  <- sort(unique(as.character(deg_stats_sex$CellType)))

deg_stats_auto <- deg_stats_auto %>%
  mutate(CellType = factor(CellType, levels = order_auto))

deg_stats_sex <- deg_stats_sex %>%
  mutate(CellType = factor(CellType, levels = order_sex))

library(ggplot2)

# Autosomal
p1 <- ggplot(deg_stats_auto, aes(x = CellType, y = prop, fill = deg_class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c( "Shared"="#377eb8", "Typical_unique"="#999999","Gonadal"="#e41a1c", "SCC"="#4daf4a"),
    breaks = c("Typical_unique", "Shared", "Gonadal", "SCC")
  ) +
  scale_y_continuous(expand = c(0, 0)) +   
  theme_bw() +
  labs(title = "Autosomal genes", y = "Proportion", fill = "Category")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 10),
    
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, "cm"),         
    #legend.key.height = unit(1, "cm"),     
    #legend.key.width  = unit(1, "cm")       
  )

# Sex-linked
p2 <- ggplot(deg_stats_sex, aes(x = CellType, y = prop, fill = deg_class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c( "Shared"="#377eb8", "Typical_unique"="#999999","Gonadal"="#e41a1c", "SCC"="#4daf4a"),
    breaks = c("Typical_unique", "Shared", "Gonadal", "SCC")
  ) +
  scale_y_continuous(expand = c(0, 0)) +   
  theme_bw() +
  labs(title = "Sex-linked genes", y = "Proportion", fill = "Category")+ # fill used for changing legend title
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), 
   
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 10),
    
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, "cm") )

ggsave("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_comparison_auto_HP.pdf", plot = p1, height = 4,width = 4)

ggsave("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_comparison_sex_HP.pdf", plot = p2, height = 4,width = 4)


deg_wide %>%
  filter(deg_class %in% c("SCC", "Gonadal")) %>%
  count(deg_class, name = "n_DEGs")
# ------------------------------------------------------------------------------

# add one without separation of autosomal and sex-linked genes

library(dplyr)

deg_stats <- deg_wide %>%
  dplyr::group_by(CellType, deg_class) %>%
  dplyr::summarise(gene_count = n(), .groups = "drop") %>%
  dplyr::group_by(CellType) %>%
  dplyr::mutate(prop = gene_count / sum(gene_count))

deg_stats <- deg_stats %>%
  mutate(deg_class = factor(deg_class, levels = c("Typical_unique", "Shared", "Gonadal", "SCC")))

order_both <- deg_stats %>%
  filter(deg_class == "Typical_unique") %>%
  arrange(prop) %>%
  pull(CellType)
# to avoid missing cell type in factors, add cell type without Normative_unique

all_celltypes <- unique(deg_stats$CellType)
order_sex_full <- c(setdiff(all_celltypes, order_both),order_both)
deg_stats <- deg_stats %>%
  mutate(CellType = factor(CellType, levels = order_sex_full))

order <- sort(unique(as.character(deg_stats$CellType)))

deg_stats <- deg_stats %>%
  mutate(CellType = factor(CellType, levels = order_auto))


p3 <- ggplot(deg_stats, aes(x = CellType, y = prop, fill = deg_class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c( "Shared"="#377eb8", "Typical_unique"="#999999","Gonadal"="#e41a1c", "SCC"="#4daf4a"),
    breaks = c("Typical_unique", "Shared", "Gonadal", "SCC")
  ) +
  scale_y_continuous(expand = c(0, 0)) +   
  theme_bw() +
  labs(title = "All genes", y = "Proportion", fill = "Category")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
   
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 10),
    
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, "cm"))

ggsave("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_comparison_wholeChrom_HP.pdf", plot = p3, height = 4,width = 4)

# p3 looks identical to autosomal genes, but the values are actually slightly different

deg_wide %>%
  filter(deg_class %in% c("SCC", "Gonadal", "Typical_unique", "Shared")) %>%
  count(deg_class, name = "n_DEGs")

# ------------------------------------------------------------------------------
# Identify top 10 up and down-regulated genes per cell type
# ------------------------------------------------------------------------------

deg_all2 <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg_all2.RDS")

deg_sig <- deg_all2 %>%
  dplyr::filter(
    !is.na(logFC),
    !is.na(adj.P.Val),
    adj.P.Val < 0.05,
    Comparison %in% c("XYT_vs_XXO", "XYT_vs_XXT", "XXT_vs_XXO")
  )

top10_up <- deg_sig %>%
  dplyr::filter(logFC > 0) %>%
  dplyr::group_by(CellType, Comparison) %>%
  dplyr::arrange(desc(logFC), .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

top10_down <- deg_sig %>%
  dplyr::filter(logFC < 0) %>%
  dplyr::group_by(CellType, Comparison) %>%
  dplyr::arrange(logFC, .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

top10_degs <- dplyr::bind_rows(
  top10_up   %>% dplyr::mutate(Direction = "Up"),
  top10_down %>% dplyr::mutate(Direction = "Down")
)

table(top10_degs$CellType,top10_degs$Comparison, top10_degs$Direction)

write.csv(
  top10_degs,
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/Top10_UpDown_DEGs_by_CellType_by_Comparison.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# make volcanoe plots
# ------------------------------------------------------------------------------


deg_all2 <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg_all2.RDS")


library(ggplot2)
library(dplyr)
library(ggrepel)

volcano_plot <- function(df, celltype, comparison,
                         fdr_cutoff = 0.05,
                         lfc_cutoff = 0.25,
                         top_n = 10) {
  
  df_sub <- df %>%
    filter(
      CellType == celltype,
      Comparison == comparison,
      !is.na(logFC),
      !is.na(adj.P.Val)
    ) %>%
    mutate(
      negLogFDR = -log10(adj.P.Val),
      Significance = case_when(
        adj.P.Val < fdr_cutoff & logFC >= lfc_cutoff  ~ "Up",
        adj.P.Val < fdr_cutoff & logFC <= -lfc_cutoff ~ "Down",
        TRUE                                   ~ "NS"
      )
    )
  
  # select top genes for labeling
  top_genes <- df_sub %>%
    filter(Significance != "NS") %>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n = top_n)
  
  ggplot(df_sub, aes(logFC, negLogFDR)) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c(
      "Up" = "#D62728",
      "Down" = "#1F77B4",
      "NS" = "grey70"
    )) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed") +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf
    ) +
    labs(
      title = paste(celltype, comparison),
      x = "log2 Fold Change",
      y = "-log10(FDR)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

celltypes <- unique(deg_all2$CellType)
comparisons <- c("XYT_vs_XXO", "XYT_vs_XXT", "XXT_vs_XXO")

pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/VolcanoPlots_All_CellTypes.pdf", width = 6, height = 6)

for (ct in celltypes) {
  for (comp in comparisons) {
    if (nrow(filter(deg_all2, CellType == ct, Comparison == comp)) > 0) {
      print(volcano_plot(deg_all2, ct, comp))
    }
  }
}

dev.off()

# ------------------------------------------------------------------------------
# plot dot plot of sex-linked gene expression changes in typical comparison
# ------------------------------------------------------------------------------

# load deg_all2

genes_of_interest <- c('Uty', 'Eif2s3y', 'Ddx3y', 'Kdm5d', 'Kdm6a', 'Kdm5c', 'Ddx3x', 'ENSRNOG00000065796', "Xist")


library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# Load Xist DEG results
# ============================================================

deg_all2Xist <- readRDS("/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Limma/deg_all2Xist.RDS")

# must have repeated all through limma analysis with Xist data after Twa et al to make deg_all2Xist.RDS
# ============================================================
# Genes of interest
# ============================================================


# ============================================================
# 1. Get all genes EXCEPT Xist from deg_all2
# ============================================================

dot_df <- deg_all2 %>%
  dplyr::filter(
    Comparison == "XYT_vs_XXO",
    gene %in% setdiff(genes_of_interest, "Xist")
  ) %>%
  dplyr::mutate(
    log10_p = -log10(adj.P.Val)
  )


# ============================================================
# 2. Get Xist from deg_all2Xist
# ============================================================

dot_Xist <- deg_all2Xist %>%
  dplyr::filter(
    Comparison == "XYT_vs_XXO",
    gene == "Xist"
  ) %>%
  dplyr::mutate(
    log10_p = -log10(adj.P.Val)
  )


# ============================================================
# 3. Combine the two DEG datasets
# ============================================================

dot_df <- dplyr::bind_rows(
  dot_df,
  dot_Xist
)


# ============================================================
# 4. Make sure every gene × cell type combination exists
# ============================================================

all_celltypes <- unique(deg_all2$CellType)

dot_df <- tidyr::expand_grid(
  gene = genes_of_interest,
  CellType = all_celltypes
) %>%
  dplyr::left_join(
    dot_df,
    by = c("gene", "CellType")
  )


# ============================================================
# 5. Order genes by regulation
# ============================================================

gene_order <- dot_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    mean_logFC = mean(logFC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(mean_logFC)) %>%
  dplyr::pull(gene)

dot_df$gene <- factor(
  dot_df$gene,
  levels = rev(gene_order)
)


# ============================================================
# 6. Dot plot
# ============================================================

p_dot <- ggplot(
  dot_df,
  aes(
    x = CellType,
    y = gene
  )
) +
  geom_point(
    aes(
      size = log10_p,
      color = logFC
    ),
    na.rm = TRUE
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  scale_size(
    range = c(1, 8)
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    x = "Cell type",
    y = "Gene",
    color = "log2FC",
    size = "-log10(adj.P.Val)"
  )


# Display
p_dot


# ============================================================
# 7. Save
# ============================================================

ggsave(
  "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Plots/DEG_dotplot_XYT_vs_XXO.pdf",
  plot = p_dot,
  width = 10,
  height = 6.5
)

################################################################################
# lookup mouse orthologs of top 10 DEGs in rat
################################################################################

# change list to look up different genes as necessary


rat_genes <- c(
  "Uty",
  "Kdm5d",
  "Ddx3",
  "Eif2s3y",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "Uty",
  "Kdm5d",
  "Ddx3",
  "Eif2s3y",
  "",
  "",
  "",
  "",
  "ENSRNOG00000065796",
  "ENSRNOG00000037911",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "ENSRNOG00000065796",
  "ENSRNOG00000037911",
  "",
  "",
  "",
  "",
  ""
)




# Your gene lookup dataframe from MS metadata
gene_lookup <- data.frame(
  Rat_gene   = rownames(MSMC),
  Mouse_gene = MSMC@assays[["RNA"]]@meta.data[["Mouse_Symbol"]],
  stringsAsFactors = FALSE
)

# To keep the original order with duplicates, do a left_join by rat gene
library(tibble)
result_ordered <- tibble(Rat_gene = rat_genes) %>%
  left_join(gene_lookup, by = "Rat_gene")

# View results with mouse gene matched, in original order

write.csv(
  result_ordered,
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/rat_to_mouse_lookup.csv",
  row.names = FALSE
)



################################################################################
# look for expression of specific genes
################################################################################

library(dplyr)
library(Seurat)

# checked Xist, Tsix - but they are not there in the dataset (Xist is after adding following Twa et al)
# same with ENSG00000229807, ENSMUSG00000086503

xist_avg <- FetchData(
MSMC,
  vars = c("Sirt1", "MyClass_Bcas1", "group")
) %>%
  group_by(MyClass_Bcas1, group) %>%
  summarise(
    mean_expr = mean(Sirt1, na.rm = TRUE),
    pct_expr  = mean(Sirt1 > 0),
    .groups = "drop"
  )

xist_avg


################################################################################
# Cell Type Proportiona Analysis
################################################################################


## ===============================
## Libraries
## ===============================
library(speckle)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

## ===============================
## 1. Metadata
## ===============================
meta <- MS@meta.data %>%
  transmute(
    celltype  = MyClass_Bcas1,
    sample_id = sample,
    group_id  = group
  )

# sanity check
table(meta$group_id, meta$sample_id)

## ===============================
## 2. Propeller analysis
## ===============================
prop <- propeller(
  clusters  = meta$celltype,
  sample    = meta$sample_id,
  group     = meta$group_id,
  transform = "logit"
)

## ===============================
## 3. Per-sample proportions (PLOTTING DATA)
## ===============================
prop_df <- meta %>%
  dplyr::count(sample_id, celltype, group_id, name = "n_cells") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(
    total_cells = sum(n_cells),
    prop = n_cells / total_cells
  ) %>%
  ungroup()

## ===============================
## 4. Propeller summary statistics
## ===============================
prop_stats <- prop %>%
  tibble::rownames_to_column("celltype") %>%
  dplyr::select(
    celltype,
    BaselineProp,
    starts_with("PropMean"),
    Fstatistic,
    P.Value,
    FDR
  )

## ===============================
## 5. Save summary table (PDF)
## ===============================
summary_grob <- tableGrob(
  prop_stats,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(cex = 0.6),
      padding = unit(c(2, 2), "mm")
    ),
    colhead = list(
      fg_params = list(cex = 0.7, fontface = "bold")
    )
  )
)

pdf(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/celltype_proportions_summary.pdf",
  width = 14,
  height = 14
)
grid.draw(summary_grob)
dev.off()

## ===============================
## 6. Significance stars
## ===============================
sig_labels <- prop %>%
  tibble::rownames_to_column("celltype") %>%
  mutate(
    sig = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    )
  ) %>%
  dplyr::select(celltype, sig)

star_pos <- prop_df %>%
  dplyr::group_by(.data$celltype) %>%
  dplyr::summarise(
    y = max(prop, na.rm = TRUE) * 1.08,
    .groups = "drop"
  )

star_df <- star_pos %>%
  dplyr::left_join(sig_labels, by = "celltype") %>%
  dplyr::filter(sig != "")



star_df <- prop_df %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    y = max(prop, na.rm = TRUE) * 1.08,
    .groups = "drop"
  ) %>%
  dplyr::left_join(sig_labels, by = "celltype") %>%
  dplyr::filter(sig != "")


## ===============================
## 7. One-page combined plot
## ===============================
group_colors <- c(
  "XXO" = "black",
  "XXT" = "black",
  "XYT" = "black"
)

p <- ggplot(
  prop_df,
  aes(x = group_id, y = prop, fill = NULL)
) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.3, alpha = 0.8) +
  facet_wrap(~ celltype, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  geom_text(
    data = star_df,
    aes(x = 2, y = y, label = sig),
    inherit.aes = FALSE,
    size = 5
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "Cell-type proportion"
  )

## ===============================
## 8. Save final figure
## ===============================
ggsave(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/celltype_proportions_all_celltypes_one_page.pdf",
  plot = p,
  width = 20,
  height = 12
)


## ===============================
## X. Group-level proportions (REQUIRED for bar plot)
## ===============================
prop_group <- prop_df %>%
  dplyr::group_by(group_id, celltype) %>%
  dplyr::summarise(
    mean_prop = mean(prop),
    .groups = "drop"
  )

## ===============================
## 10. Stacked bar plot
## ===============================
p_bar <- ggplot(prop_group, aes(x = group_id, y = mean_prop, fill = celltype)) +
  geom_bar(
    stat = "identity",
    position = "fill",
    width = 0.7,
    color = "black",     # black borders between segments
    linewidth = 0.5      # adjust border thickness
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Cell-type composition (%)",
    fill = "Cell type"
  )

## ===============================
## 11. Save bar plot
## ===============================
ggsave(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/celltype_composition_barplot.pdf",
  plot = p_bar,
  width = 10,
  height = 6
)

# by sample

library(dplyr)

prop_sample <- MS@meta.data %>%
  group_by(sample, group, MyClass_Bcas1) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample, group) %>%
  mutate(mean_prop = n / sum(n)) %>%
  ungroup()

library(ggplot2)

library(dplyr)

prop_sample <- prop_sample %>%
  mutate(
    sample_id = sub(".*(MS[0-9]+).*", "\\1", sample),
    sample_label = paste(sample_id, group, sep = "\n")
  )

p_bar <- ggplot(prop_sample,
                aes(x = sample_label, y = mean_prop, fill = MyClass_Bcas1)) +
  geom_bar(
    stat = "identity",
    position = "fill",
    width = 0.7,
    color = "black",
    linewidth = 0.5
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 25)
  ) +
  labs(
    x = NULL,
    y = "Cell-type composition (%)",
    fill = "Cell type"
  )

ggsave(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/celltype_composition_barplot_by_sample.pdf",
  plot = p_bar,
  width = 12,
  height = 8
)

################################################################################
# clustering DEGs by correlation analysis
################################################################################

deg_all2 <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg_all2.RDS")

# ==============================================================================
# Libraries
# ==============================================================================
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ==============================================================================
# 1. Rename & relabel comparison → Effect
# ==============================================================================
deg <- deg_all2 %>%
  dplyr::rename(Effect = Comparison) %>%
  dplyr::mutate(
    Effect = recode(
      Effect,
      "XYT_vs_XXO" = "Typical",
      "XXT_vs_XXO" = "Gonadal",
      "XYT_vs_XXT" = "SCC"
    )
  )

# ==============================================================================
# 2. Aggregate logFC per gene × CellType × Effect
# ==============================================================================
deg_unique <- deg %>%
  dplyr::group_by(gene, CellType, Effect) %>%
  dplyr::summarise(
    logFC = mean(logFC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Group = paste(CellType, Effect, sep = " | ")
  )

# ==============================================================================
# 3. Wide matrix: genes × (CellType | Effect)
# ==============================================================================
logFC_mat <- deg_unique %>%
  dplyr::select(gene, Group, logFC) %>%
  pivot_wider(names_from = Group, values_from = logFC)

logFC_mat <- as.data.frame(logFC_mat)
rownames(logFC_mat) <- logFC_mat$gene
logFC_mat$gene <- NULL
logFC_mat[] <- lapply(logFC_mat, as.numeric)
logFC_mat <- as.matrix(logFC_mat)

# ==============================================================================
# 4. Correlation matrix (shared genes only)
# ==============================================================================
cor_mat <- cor(
  logFC_mat,
  use = "pairwise.complete.obs",
  method = "pearson"
)

# ==============================================================================
# 5. Correlation p-value matrix
# ==============================================================================
get_cor_pmat <- function(mat) {
  n <- ncol(mat)
  pmat <- matrix(NA, n, n)
  colnames(pmat) <- rownames(pmat) <- colnames(mat)
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) {
        pmat[i, j] <- 1
      } else {
        pmat[i, j] <- cor.test(mat[, i], mat[, j])$p.value
      }
    }
  }
  pmat
}

p_mat <- get_cor_pmat(logFC_mat)

# ==============================================================================
# 6. Significance stars
# ==============================================================================
get_star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

ann_star <- get_star(p_mat)

# ==============================================================================
# 7. Annotation data (CellType + Effect)
# ==============================================================================
annotation_df <- data.frame(
  Group = colnames(cor_mat)
) %>%
  separate(Group, into = c("CellType", "Effect"), sep = " \\| ")

rownames(annotation_df) <- colnames(cor_mat)

# ==============================================================================
# 8. Annotation colors
# ==============================================================================
ann_colors <- list(
  CellType = c(
    "Mature Oligodendrocytes" = "darkred",
    "Immature Oligodendrocytes" = "red",
    "GABAergic Neurons (1)" = "purple",
    "GABAergic Neurons (2)" = "lavender",
    "Astrocytes" = "yellow",
    "OPCs" = "orange",
    "Glutamatergic Neurons" = "green",
    "Cholinergic Neurons" = "blue"
  ),
  Effect = c(
    "Typical"  = "black",
    "Gonadal"  = "darkgray",
    "SCC"      = "lightgray"
  )
)


# ==============================================================================
# 9. ComplexHeatmap
# ==============================================================================
ht <- Heatmap(
  cor_mat,
  name = "Pearson r",
  
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  
  # -----------------------
  # Dendrogram size (KEY FIX)
  # -----------------------
  column_dend_height = unit(3, "cm"),
  row_dend_width     = unit(3, "cm"),
  
  # Optional: thicker dendrogram lines
  column_dend_gp = gpar(lwd = 3),
  row_dend_gp    = gpar(lwd = 3),
  
  # -----------------------
  # Axis labels
  # -----------------------
  row_names_gp    = gpar(fontsize = 40),
  column_names_gp = gpar(fontsize = 40),
  
  # -----------------------
  # Title
  # -----------------------
  column_title = "logFC Correlations by Cell Type and Effect\n(Autosomal + Sex-linked)",
  column_title_gp = gpar(fontsize = 40, fontface = "bold"),
  
  # -----------------------
  # Heatmap legend
  # -----------------------
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp  = gpar(fontsize = 36, fontface = "bold"),
    labels_gp = gpar(fontsize = 36),
    at = c(-1, 0, 1)
  ),
  
  # -----------------------
  # Top annotation
  # -----------------------
  top_annotation = HeatmapAnnotation(
    df = annotation_df,
    col = ann_colors,
    simple_anno_size = unit(1.2, "cm"),
    annotation_name_gp = gpar(fontsize = 40, fontface = "bold"),
    annotation_legend_param = list(
      title_gp  = gpar(fontsize = 38, fontface = "bold"),
      labels_gp = gpar(fontsize = 36)
    )
  ),
  
  # -----------------------
  # Left annotation
  # -----------------------
  left_annotation = rowAnnotation(
    df = annotation_df,
    col = ann_colors,
    simple_anno_size = unit(1.2, "cm"),
    annotation_name_gp = gpar(fontsize = 40, fontface = "bold"),
    annotation_legend_param = list(
      title_gp  = gpar(fontsize = 38, fontface = "bold"),
      labels_gp = gpar(fontsize = 36)
    )
  ),
  
  show_row_names    = TRUE,
  show_column_names = FALSE,
  
  # -----------------------
  # Cell text (stars)
  # -----------------------
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      ann_star[i, j],
      x, y,
      gp = gpar(fontsize = 20, fontface = "bold")
    )
  }
)


# ==============================================================================
# 10. Save
# ==============================================================================
pdf(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/Correlation_CellType_Effect.pdf",
  width = 42, height = 26
)
draw(
  ht,
  annotation_legend_side = "bottom",
  heatmap_legend_side    = "bottom",
  merge_legend           = TRUE,
  padding = unit(c(20, 5, 10, 250), "mm"))
dev.off()

# cluster groups in heatmap and name them

deg_all2 <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg_all2.RDS")

C1 <- c("Mature Oligodendrocytes | Typical", "Immature Oligodendrocytes | Typical", 
        "Mature Oligodendrocytes | SCC", "Immature Oligodendrocytes | SCC")
# Oligodendrocytes | Typical & SCC


C2 <- c("Astrocytes | Typical", "Astrocytes | SCC")
# Astrocytes | Typical & SCC


C3 <- c("Glutamatergic Neurons | Typical", "Glutamatergic Neurons | SCC")
# Glutamatergic Neurons | Typical & SCC

C4 <- c("Mature Oligodendrocytes | Gonadal", "Immature Oligodendrocytes | Gonadal")
# Oligodendrocytes | Gonadal

C5 <- c("OPCs | Typical", "OPCs | Gonadal")
# OPCs | Typical & Gonadal

C6 <- "Astrocytes | Gonadal"
# Astrocytes | Gonadal

C7 <- c("GABAergic Neurons (1) | Typical", "GABAergic Neurons (1) | SCC")
# GABAergic Neurons (1) | Typical & SCC


C8 <- "OPCs | SCC"
# OPCs | SCC

C9 <- c("Cholinergic Neurons | Typical", "Cholinergic Neurons | SCC")
# Cholingeric Neurons | Typical & SCC

C10 <- c("GABAergic Neurons (1) | Gonadal", "Glutamatergic Neurons | Gonadal", 
         "Cholinergic Neurons | Gonadal")
# Neurons | Gonadal

C11 <- c("GABAergic Neurons (2) | Typical", "GABAergic Neurons (2) | SCC")
# GABAergic Neurons (2) | Typical & SCC

C12 <- "GABAergic Neurons (2) | Gonadal"
# GABAergic Neurons (2) | Gonadal

clusters <- list(
  C1 = C1, C2 = C2, C3 = C3, C4 = C4,
  C5 = C5, C6 = C6, C7 = C7, C8 = C8, C9 = C9, C10 = C10, C11 = C11, C12 = C12
)


# ------------------------------------------------------------------------------
# Clusterprofiler
# ------------------------------------------------------------------------------


library(dplyr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(purrr)

deg <- deg_all2 %>%
  dplyr::rename(Effect = Comparison) %>%
  dplyr::mutate(
    Effect = recode(
      Effect,
      "XYT_vs_XXO" = "Typical",
      "XXT_vs_XXO" = "Gonadal",
      "XYT_vs_XXT" = "SCC"
    ),
    Group = paste(CellType, Effect, sep = " | ")
  )

assign_cluster <- function(group, clusters) {
  cl <- names(clusters)[purrr::map_lgl(clusters, ~ group %in% .x)]
  if (length(cl) == 0) NA_character_ else cl
}

deg <- deg %>%
  mutate(
    Cluster = purrr::map_chr(Group, assign_cluster, clusters = clusters)
  ) %>%
  filter(!is.na(Cluster))

unique(deg$Cluster)


# save deg ---------------------------------------------------------------------
# saveRDS(deg, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg.RDS")
# deg <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg.RDS")
# ------------------------------------------------------------------------------

################################################################################
# ORA analysis
################################################################################

library(clusterProfiler)
library(org.Rn.eg.db)
library(AnnotationDbi)

run_ora_cluster <- function(df, cluster_name) {
  
  message("Running ORA (BP): ", cluster_name)
  
  sub <- df %>%
    filter(
      Cluster == cluster_name,
      adj.P.Val < 0.05
    )
  
  if (nrow(sub) == 0) return(NULL)
  
  # ---------------------------
  # 1. Split gene IDs
  # ---------------------------
  genes_symbol <- sub$gene[!grepl("^ENS", sub$gene)]
  genes_ens    <- sub$gene[ grepl("^ENS", sub$gene)]
  
  # ---------------------------
  # 2. SYMBOL → ENTREZ
  # ---------------------------
  map_sym <- suppressMessages(
    bitr(unique(genes_symbol),
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = org.Rn.eg.db)
  )
  
  # ---------------------------
  # 3. ENSEMBL → ENTREZ
  # ---------------------------
  valid_ens <- intersect(
    unique(genes_ens),
    keys(org.Rn.eg.db, keytype = "ENSEMBL")
  )
  
  map_ens <- suppressMessages(
    bitr(valid_ens,
         fromType = "ENSEMBL",
         toType   = "ENTREZID",
         OrgDb    = org.Rn.eg.db)
  )
  
  genes <- unique(c(map_sym$ENTREZID, map_ens$ENTREZID))
  
  if (length(genes) == 0) return(NULL)
  
  # ---------------------------
  # 4. Universe (all tested genes)
  # ---------------------------
  all_genes <- unique(df$gene)
  
  uni_sym <- all_genes[!grepl("^ENS", all_genes)]
  uni_ens <- all_genes[ grepl("^ENS", all_genes)]
  
  uni_map_sym <- suppressMessages(
    bitr(uni_sym,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = org.Rn.eg.db)
  )
  
  valid_uni_ens <- intersect(
    uni_ens,
    keys(org.Rn.eg.db, keytype = "ENSEMBL")
  )
  
  uni_map_ens <- suppressMessages(
    bitr(valid_uni_ens,
         fromType = "ENSEMBL",
         toType   = "ENTREZID",
         OrgDb    = org.Rn.eg.db)
  )
  
  universe <- unique(c(
    uni_map_sym$ENTREZID,
    uni_map_ens$ENTREZID
  ))
  
  # ---------------------------
  # 5. ORA (BP)
  # ---------------------------
  suppressMessages(
    enrichGO(
      gene          = genes,
      universe      = universe,
      OrgDb         = org.Rn.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      readable      = TRUE
    )
  )
}

library(purrr)

ego_results <- tibble(
  Cluster = names(clusters)
) %>%
  mutate(
    ego = map(Cluster, ~ run_ora_cluster(deg, .x))
  )

library(tidyr)

ego_table <- ego_results %>%
  filter(!map_lgl(ego, is.null)) %>%
  mutate(result = map(ego, as.data.frame)) %>%
  unnest(result)

ego_table_clean <- ego_table %>%
  dplyr::select(-ego)

write.csv(
  ego_table_clean,
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/GO_BP_ORA_all_clusters.csv",
  row.names = FALSE
)

library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(tidyr)

plot_one_cluster <- function(ego, cluster) {
  
  if (is.null(ego)) return(NULL)
  
  df <- as.data.frame(ego)
  if (nrow(df) == 0) return(NULL)
  
  dotplot(ego, showCategory = 15) +
    ggtitle(paste(cluster, "GO BP")) +
    theme(plot.title = element_text(face = "bold"))
}

plots <- ego_results %>%
  mutate(
    plot = map2(
      ego,
      Cluster,
      plot_one_cluster
    )
  ) %>%
  filter(!map_lgl(plot, is.null))

plots %>%
  pwalk(function(Cluster, plot, ...) {
    
    fname <- paste0(
      "GO_BP__",
      gsub("[^A-Za-z0-9]", "_", Cluster),
      ".png"
    )
    
    ggsave(
      filename = file.path(
        "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/",
        fname
      ),
      plot  = plot,
      width = 8,
      height = 6,
      dpi   = 300
    )
  })

cluster_labels <- c(
  "C1" = "Oligodendrocytes | Typical & SCC",
  "C2" = "Astrocytes | Typical & SCC",
  "C3" = "Glutamatergic Neurons | Typical & SCC",
  "C4" = "Oligodendrocytes | Gonadal",
  "C5" = "OPCs | Typical & Gonadal",
  "C6" = "Astrocytes | Gonadal",
  "C7" = "GABAergic Neurons (1) | Typical & SCC",
  "C8" = "OPCs | SCC", 
  "C9" = "Cholinergic Neurons | Typical & SCC",
  "C10" = "Neurons | Gonadal",
  "C11" = "GABAergic Neurons (2) | Typical & SCC",
  "C12" = "GABAergic Neurons (2) | Gonadal"
)

ego_df <- ego_results %>%
  filter(!map_lgl(ego, is.null)) %>%
  transmute(
    Cluster,
    result = map(ego, as.data.frame)
  ) %>%
  unnest(result)

top_terms <- ego_df %>%
  dplyr::group_by(Description) %>%
  dplyr::summarise(
    min_p = min(`p.adjust`),
    .groups = "drop"
  ) %>%
  dplyr::arrange(min_p) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(Description)

ego_df <- ego_df %>%
  filter(Description %in% top_terms)

ego_df <- ego_df %>%
  dplyr::mutate(
    Cluster = dplyr::recode(Cluster, !!!cluster_labels)
  )

combined_dotplot <- ggplot(
  ego_df,
  aes(
    x     = Cluster,
    y     = Description,
    size  = -log10(p.adjust),
    color = Count
  )
) +
  geom_point() +
  scale_color_continuous(
    low  = "red",
    high = "blue",
    name = "Gene count"
  ) +
  scale_size(range = c(2, 8), name = "-log10(adj.P)") +
  labs(
    x = NULL,
    y = "GO Biological Process",
    title = "GO BP enrichment across clusters"
  ) +
  theme_bw() +
  theme(
    axis.text.y   = element_text(size = 9),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "grey90"),
    plot.title    = element_text(face = "bold")
  )

ggsave(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/GO_BP_COMBINED_DOTPLOT_CLUSTERS.png",
  combined_dotplot,
  width  = 14,
  height = 8,
  dpi    = 300
)


################################################################################
# rrvgo reduction of ORA results
################################################################################

library(rrvgo)
library(clusterProfiler)
library(GOSemSim)
library(org.Rn.eg.db)
library(dplyr)
library(purrr)

reduce_go_terms <- function(ego,
                            ont = "BP",
                            threshold = 0.7) {
  
  if (is.null(ego)) return(NULL)
  
  df <- as.data.frame(ego)
  if (nrow(df) < 2) return(df)
  
  ## Semantic data
  semData <- godata(
    OrgDb   = org.Rn.eg.db,
    ont     = ont,
    keytype = "ENTREZID"
  )
  
  ## Similarity matrix
  simMatrix <- mgoSim(
    df$ID,
    df$ID,
    semData,
    measure = "Wang",
    combine = NULL
  )
  
  ## Scores
  scores <- -log10(df$p.adjust)
  names(scores) <- df$ID
  
  ## Reduce
  reduced <- reduceSimMatrix(
    simMatrix,
    scores    = scores,
    threshold = threshold,
    orgdb     = org.Rn.eg.db,
    keytype   = "ENTREZID",
    children  = TRUE
  )
  
  ## ---- IMPORTANT PART ----
  ## The representative GO term is always the FIRST column
  rep_col <- colnames(reduced)[1]
  
  reduced %>%
    dplyr::rename(GO_ID = all_of(rep_col)) %>%
    left_join(
      df,
      by = c("GO_ID" = "ID")
    )
}

ego_reduced <- ego_results %>%
  mutate(
    reduced = map(ego, reduce_go_terms)
  )

ego_reduced_table <- ego_reduced %>%
  filter(!map_lgl(reduced, is.null)) %>%
  dplyr::select(Cluster, reduced) %>%
  unnest(reduced) %>%
  distinct(Cluster, GO_ID, .keep_all = TRUE) %>%
  dplyr::select(
    Cluster,
    GO_ID,
    Description,
    p.adjust,
    Count,
    GeneRatio
  )


write.csv(
  ego_reduced_table,
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/GO_BP_ORA_REDUCED_CLUSTERS.csv",
  row.names = FALSE
)

# plot rrvgo results ------

library(dplyr)
library(ggplot2)
library(forcats)
library(purrr)

top_terms <- ego_reduced_table %>%
  dplyr::group_by(Description) %>%
  dplyr::summarise(min_p = min(p.adjust), .groups = "drop") %>%
  dplyr::arrange(min_p) %>%
  dplyr::slice_head(n = 25) %>%
  dplyr::pull(Description)

plot_df <- ego_reduced_table %>%
  
  filter(
    p.adjust < 0.05,
    Description %in% top_terms
  ) %>%
  
  mutate(
    neglogP = -log10(p.adjust)
  ) %>%
  
  group_by(Description) %>%
  mutate(
    n_clusters = n_distinct(Cluster)
  ) %>%
  ungroup() %>%
  
  arrange(n_clusters, Description) %>%
  
  mutate(
    Description = fct_reorder(Description, n_clusters, .fun = max)
  )


rrvgo_combined_dotplot <- ggplot(
  plot_df,
  aes(
    x     = Cluster,
    y     = Description,
    size  = neglogP,
    color = Count
  )
) +
  geom_point(alpha = 0.85) +
  
  scale_color_continuous(
    low  = "red",
    high = "blue",
    name = "DEG count"
  ) +
  
  scale_size(
    range = c(2, 8),
    name  = "-log10(adj.P)"
  ) +
  
  labs(
    x = NULL,
    y = "Representative GO Biological Process",
    title = "Reduced GO BP enrichment across clusters"
  ) +
  
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "grey90"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/GO_BP_RRVGO_COMBINED_DOTPLOT_CLUSTERS.png",
  plot     = rrvgo_combined_dotplot,
  width    = 14,
  height   = 8,
  dpi      = 300
)

# ------------------------------------------------------------------------------
# add logFC data after ORA and rrvgo
# ------------------------------------------------------------------------------ 

deg_fc <- deg %>%
  filter(!is.na(logFC)) %>%
  group_by(Cluster, gene) %>%
  summarise(
    logFC = mean(logFC),
    .groups = "drop"
  )

library(tidyr)

deg_fc <- deg %>%
  dplyr::select(Cluster, gene, logFC) %>%   # force columns
  dplyr::filter(!is.na(logFC)) %>%
  dplyr::group_by(Cluster, gene) %>%
  dplyr::summarise(
    logFC = mean(logFC),
    .groups = "drop"
  )
library(stringr)
library(tidyr)
library(dplyr)

ora_long <- ego_results %>%
  filter(!map_lgl(ego, is.null)) %>%
  transmute(
    Cluster,
    result = map(ego, as.data.frame)
  ) %>%
  unnest(result) %>%
  dplyr::select(
    Cluster,
    ID,
    Description,
    geneID
  ) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(gene = geneID)


go_mean_logFC <- ora_long %>%
  dplyr::inner_join(
    deg_fc,
    by = c("Cluster", "gene")
  ) %>%
  dplyr::group_by(Cluster, ID, Description) %>%
  dplyr::summarise(
    mean_logFC = mean(logFC, na.rm = TRUE),
    n_genes    = dplyr::n(),
    .groups    = "drop"
  )


go_mean_logFC <- ora_long %>%
  dplyr::inner_join(
    deg_fc,
    by = c("Cluster", "gene")
  ) %>%
  dplyr::group_by(Cluster, ID, Description) %>%
  dplyr::summarise(
    mean_logFC = mean(logFC),
    n_genes    = n(),
    .groups    = "drop"
  )



ego_reduced_table <- ego_reduced_table %>%
  left_join(
    go_mean_logFC,
    by = c(
      "Cluster",
      "GO_ID" = "ID"
    )
  )

colnames(ego_reduced_table)

ego_reduced_table <- ego_reduced_table %>%
  dplyr::rename(
    Description = Description.x
  ) %>%
  dplyr::select(-Description.y)

colnames(ego_reduced_table)

ego_reduced_table %>%
  dplyr::select(Cluster, Description, mean_logFC, Count, n_genes) %>%
  head(10)

# --- plotting

library(dplyr)

top_terms <- ego_reduced_table %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::group_by(Description) %>%
  dplyr::summarise(
    min_p = min(p.adjust),
    .groups = "drop"
  ) %>%
  dplyr::arrange(min_p) %>%
  dplyr::slice_head(n = 25) %>%
  dplyr::pull(Description)

library(forcats)

plot_df <- ego_reduced_table %>%
  
  filter(
    p.adjust < 0.05,
    Description %in% top_terms
  ) %>%
  
  mutate(
    neglogP = -log10(p.adjust)
  ) %>%
  
  group_by(Description) %>%
  mutate(
    n_clusters = n_distinct(Cluster)
  ) %>%
  ungroup() %>%
  
  arrange(n_clusters, Description) %>%
  
  mutate(
    Description = fct_reorder(Description, n_clusters, .fun = max)
  )

library(ggplot2)

plot_df <- plot_df %>%
  dplyr::mutate(
    Cluster = dplyr::recode(Cluster, !!!cluster_labels)
  )

rrvgo_dotplot <- ggplot(
  plot_df,
  aes(
    x     = Cluster,
    y     = Description,
    size  = neglogP,
    color = mean_logFC
  )
) +
  geom_point(alpha = 0.85) +
  
  scale_color_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    name     = "Mean logFC"
  ) +
  
  scale_size(
    range = c(2.5, 9),
    name  = expression(-log[10]~"(adj. P)")
  ) +
  
  labs(
    x = NULL,
    y = "Representative GO Biological Process",
    title = "Top 25 reduced GO BP terms across DEG clusters",
    subtitle = "ORA + rrvgo; dot color reflects mean logFC"
  ) +
  
  theme_bw() +
  theme(
    axis.text.y        = element_text(size = 16),
    axis.text.x        = element_text(size = 16, angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "grey90"),
    plot.title         = element_text(size = 16, face = "bold"),
    legend.title       = element_text(size = 16, face = "bold")
  )
ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/GO_BP_RRVGO_TOP25_GLOBAL.png",
  plot     = rrvgo_dotplot,
  width    = 16,
  height   = 12,
  dpi      = 300
)


################################################################################
# Enrichr
################################################################################

deg <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg.RDS")

library(tidyverse)
library(readr)

convert <- read_csv(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/HP/MapMyCells/mammalian_orthologs_20231113.csv",
  show_col_types = FALSE
)

colnames(convert)

map_rat_symbol <- convert %>%
  filter(
    !is.na(Rat_Symbol),
    !is.na(Mouse_Symbol)
  ) %>%
  distinct(Rat_Symbol, .keep_all = TRUE) %>%
  dplyr::select(
    gene = Rat_Symbol,
    mouse_gene = Mouse_Symbol
  )

map_rat_ens <- convert %>%
  filter(
    !is.na(Rat_EnsemblID),
    !is.na(Mouse_Symbol)
  ) %>%
  distinct(Rat_EnsemblID, .keep_all = TRUE) %>%
  dplyr::select(
    gene = Rat_EnsemblID,
    mouse_gene = Mouse_Symbol
  )

deg_mouse <- deg %>%
  
  # Try symbol-based mapping first
  left_join(map_rat_symbol, by = "gene") %>%
  
  # Try Ensembl-based mapping for remaining
  left_join(
    map_rat_ens,
    by = "gene",
    suffix = c("", "_ens")
  ) %>%
  
  # Final mouse gene symbol
  mutate(
    gene_mouse = coalesce(mouse_gene, mouse_gene_ens)
  ) %>%
  
  dplyr::select(-mouse_gene, -mouse_gene_ens)

deg_clean <- deg_mouse %>%
  
  filter(
    adj.P.Val < 0.05,
    !is.na(gene_mouse),
    !is.na(Cluster)
  ) %>%
  
  group_by(Cluster, gene_mouse) %>%
  
  slice_max(
    order_by = abs(logFC),
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup()

deg_clean %>%
  summarise(
    total_rows = n(),
    n_clusters = n_distinct(Cluster),
    n_genes    = n_distinct(gene_mouse)
  )

deg_clean %>%
  count(Cluster, sort = TRUE)


out_dir <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Enrichr/Enrichr_inputs"

deg_clean %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(genes = list(unique(gene_mouse)), .groups = "drop") %>%
  pwalk(function(Cluster, genes) {
    writeLines(
      genes,
      file.path(out_dir, paste0("genes_", Cluster, ".txt"))
    )
  })

# manually input genes from each Enrichr Input file in output directory into Enrichr website
# https://maayanlab.cloud/enrichr/

# --- read in results from EnrichR


library(tidyverse)

enrichr_dir <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Enrichr/Enrichr_results"

jaspar_files <- list.files(
  enrichr_dir,
  pattern = "^C[1-8]\\.txt$",
  full.names = TRUE
)

jaspar_all <- map_dfr(jaspar_files, function(f) {
  
  cluster <- tools::file_path_sans_ext(basename(f))  # "C1", "C2", ...
  
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(Cluster = cluster)
})

glimpse(jaspar_all)

jaspar_clean <- jaspar_all %>%
  dplyr::rename(
    TF        = Term,
    pvalue    = `P-value`,
    p.adjust  = `Adjusted P-value`,
    oddsratio = `Odds Ratio`,
    combined  = `Combined Score`,
    genes     = Genes
  ) %>%
  dplyr::filter(!is.na(p.adjust))

jaspar_clean <- jaspar_clean %>%
  mutate(
    # Number of overlapping genes
    n_genes = as.numeric(str_extract(Overlap, "^[0-9]+")),
    
    # Clean TF names (optional but recommended)
    TF = str_remove(TF, "_.*")
  )

jaspar_sig <- jaspar_clean %>%
  filter(p.adjust < 0.05)

# Significant TFs per cluster
jaspar_sig %>%
  dplyr::count(Cluster, sort = TRUE)

# Most significant TFs globally
jaspar_sig %>%
  group_by(TF) %>%
  summarise(min_p = min(p.adjust)) %>%
  arrange(min_p) %>%
  head(10)

write_csv(
  jaspar_sig,
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Enrichr/output.csv"
)

# retrieve logFC data

library(tidyr)
library(dplyr)
library(stringr)

jaspar_long <- jaspar_clean %>%
  separate_rows(genes, sep = ";") %>%
  mutate(
    genes = str_trim(genes)
  )

jaspar_fc <- jaspar_long %>%
  left_join(
    deg_clean %>%
      dplyr::select(Cluster, gene_mouse, logFC),
    by = c("Cluster", "genes" = "gene_mouse")
  )

tf_logfc <- jaspar_fc %>%
  dplyr::filter(!is.na(logFC)) %>%
  dplyr::group_by(Cluster, TF) %>%
  dplyr::summarise(
    mean_logFC = mean(logFC),
    n_genes_fc = n(),
    .groups = "drop"
  )

jaspar_plot_df <- jaspar_clean %>%
  left_join(
    tf_logfc,
    by = c("Cluster", "TF")
  ) %>%
  filter(!is.na(mean_logFC))

# --- plot

top_tfs <- jaspar_plot_df %>%
  dplyr::group_by(TF) %>%
  dplyr::summarise(min_p = min(p.adjust), .groups = "drop") %>%
  dplyr::arrange(min_p) %>%
  dplyr::slice_head(n = 25) %>%
  dplyr::pull(TF)

jaspar_plot_df <- jaspar_plot_df %>%
  filter(TF %in% top_tfs)

jaspar_plot_df <- jaspar_plot_df %>%
  group_by(TF) %>%
  mutate(n_clusters = n_distinct(Cluster)) %>%
  ungroup() %>%
  arrange(n_clusters) %>%
  mutate(
    TF = forcats::fct_reorder(TF, n_clusters, .fun = max)
  )

library(ggplot2)

jaspar_plot_df <- jaspar_plot_df %>%
  dplyr::mutate(
    Cluster = dplyr::recode(Cluster, !!!cluster_labels)
  )

tf_dotplot <- ggplot(
  jaspar_plot_df,
  aes(
    x     = Cluster,
    y     = TF,
    size  = -log10(p.adjust),
    color = mean_logFC
  )
) +
  geom_point(alpha = 0.85) +
  
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Mean logFC"
  ) +
  
  scale_size(
    range = c(2, 8),
    name  = "-log10(adj.P)"
  ) +
  
  labs(
    x = NULL,
    y = "Transcription Factor (JASPAR)",
    title = "JASPAR TF enrichment across clusters\nColored by mean logFC of target genes"
  ) +
  
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "grey90"),
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/JASPAR_TF_DOTPLOT_mean_logFC.png",
  plot     = tf_dotplot,
  width    = 6,
  height   = 12,
  dpi      = 300
)



################################################################################
# SCING
################################################################################

# install SCING

# bash

# git clone https://github.com/XiaYangLabOrg/SCING.git
# cd SCING

# conda env create -n scing --file install/scing.environment.yml  
# conda activate scing
# pip install pyitlib  
# pip install -e .

# retrieve metacell object data and prepare for SCING in python by converting to H5AD

MSMC <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MSMC.RDS")


library(Seurat)
library(Matrix)

DefaultAssay(MSMC) <- "RNA"

counts <- GetAssayData(MSMC, layer = "counts")
data   <- GetAssayData(MSMC, layer = "data")

meta <- MSMC@meta.data

MSMC_clean <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  assay = "RNA"
)

# Put log-normalized data back
MSMC_clean <- SetAssayData(
  MSMC_clean,
  layer = "data",
  new.data = data
)

DefaultAssay(MSMC_clean) <- "RNA"

validObject(MSMC_clean)  # MUST PASS

library(MuDataSeurat)

out_dir <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/adata_by_celltype"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

celltypes <- unique(MSMC_clean@meta.data$MyClass_Bcas1)

for (ct in celltypes) {
  
  message("Processing ", ct)
  
  sub <- subset(
    MSMC_clean,
    subset = MyClass_Bcas1 == ct
  )
  
  # Skip if no cells
  if (ncol(sub) == 0) {
    message("  -> skipped (no cells)")
    next
  }
  
  # Slim object
  sub <- DietSeurat(
    sub,
    assays = "RNA",
    layers = c("counts", "data"),
    features = rownames(sub)
  )
  
  validObject(sub)
  
  MuDataSeurat::WriteH5AD(
    sub,
    file = file.path(
      out_dir,
      paste0("MSMC_", ct, ".h5ad")
    ),
    assay = "RNA"
  )
}

# create .sh and .py scripts to run SCING

# /u/scratch/v/vturnbil/GSU_FCG/Restart/MS/SCING/qsub_SCING_MSMC.sh

# /u/scratch/v/vturnbil/GSU_FCG/Restart/MS/SCING/run_scing_MSMC.py




# will need /u/scratch/v/vturnbil/GSU_FCG/Restart/MS/SCING/filelist.txt




# then run .sh and .R scripts to run wKDA analysis after SCING

# /u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/qsub_wKDA_MSMC.sh

# /u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/R_Files/wKDA_after_SCING.R



################################################################################
# run MDF in python to start Mergeomics analysis
################################################################################

# run shell_filestxt_gwas.sh (and filestxt_gwas.py) to make newfiles.txt with list
# of GWAS dataset tsv files and columns to read from. Manually add skipped files
# (look at log to check column names and see skipped files). Then rename newfiles.txt
# as files.txt

# run shell_preprocess_gwas.sh (and preprocess_gwas.py) to preprocess GWAS tsv files into txt

# run MDF.bash (with mdprune) to complete MDF and generate genes.txt and marker.txt files



# convert rat genes to human orthologs and then make mod and inf files

deg <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg_all2.RDS")

library(tidyverse)
library(readr)

convert <- read_csv(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/HP/MapMyCells/mammalian_orthologs_20231113.csv",
  show_col_types = FALSE
)

colnames(convert)

# Mouse → Human mapping (symbol-based)
map_rat_symbol <- convert %>%
  filter(
    !is.na(Rat_Symbol),
    !is.na(Human_Symbol)
  ) %>%
  distinct(Rat_Symbol, .keep_all = TRUE) %>%
  dplyr::select(
    gene = Rat_Symbol,
    human_gene = Human_Symbol
  )

map_rat_ens <- convert %>%
  filter(
    !is.na(Rat_EnsemblID),
    !is.na(Human_Symbol)
  ) %>%
  distinct(Rat_EnsemblID, .keep_all = TRUE) %>%
  dplyr::select(
    gene = Rat_EnsemblID,
    human_gene = Human_Symbol
  )

deg_human <- deg %>%
  
  # Try symbol-based mapping first
  left_join(map_rat_symbol, by = "gene") %>%
  
  # Try Ensembl-based mapping for remaining
  left_join(
    map_rat_ens,
    by = "gene",
    suffix = c("", "_ens")
  ) %>%
  
  # Final mouse gene symbol
  mutate(
    gene_human = coalesce(human_gene, human_gene_ens)
  ) %>%
  
  dplyr::select(-human_gene, -human_gene_ens)

deg_clean <- deg_human %>%
  
  filter(
    adj.P.Val < 0.05,
    !is.na(gene_human),
    !is.na(CellType)
  ) %>%
  
  group_by(CellType, gene_human) %>%
  
  slice_max(
    order_by = abs(logFC),
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup()

library(dplyr)

# -------------------------------
# 1. CLEAN INPUT
# -------------------------------

deg_clean <- deg_clean %>%
  filter(!is.na(gene_human), gene_human != "")

# -------------------------------
# 2. SELECT TOP GENES
# per CellType × Comparison
# -------------------------------

deg_clean <- deg_clean %>%
  arrange(adj.P.Val) %>%
  group_by(CellType, Comparison) %>%
  slice_head(n = 500) %>%
  ungroup()

# -------------------------------
# 3. SUMMARY CHECK
# -------------------------------

deg_clean %>%
  summarise(
    total_rows = n(),
    n_celltypes = n_distinct(CellType),
    n_comparisons = n_distinct(Comparison),
    n_genes = n_distinct(gene_human)
  )

# -------------------------------
# 4. BUILD MODULE FILE
# modules = CellType × Comparison
# -------------------------------

mod_df <- deg_clean %>%
  mutate(module = paste(CellType, Comparison, sep = "_")) %>%
  dplyr::select(module, gene = gene_human) %>%
  distinct()

# -------------------------------
# 5. CLEAN MODULE NAMES
# -------------------------------

mod_df <- mod_df %>%
  mutate(
    module = module %>%
      gsub(" ", "_", .) %>%
      gsub("[^A-Za-z0-9_]", "", .)
  )

colnames(mod_df) <- c("MODULE","GENE")

# -------------------------------
# 6. FILTER SMALL MODULES
# -------------------------------

mod_df <- mod_df %>%
  group_by(MODULE) %>%
  filter(n() >= 10) %>%
  ungroup()

# -------------------------------
# 7. WRITE MODULE FILE
# -------------------------------

write.table(
  mod_df,
  file = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/ratDE_human.mod",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# -------------------------------
# 8. BUILD INFO FILE
# -------------------------------

inf_df <- mod_df %>%
  distinct(MODULE) %>%
  mutate(
    source = "DEG",
    description = "MS"
  ) %>%
  dplyr::select(MODULE, source, description)

colnames(inf_df) <- c("MODULE","SOURCE","DESCR")

# -------------------------------
# 9. WRITE INFO FILE
# -------------------------------

write.table(
  inf_df,
  file = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/ratDE_human.inf",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# -------------------------------
# 10. FINAL CHECK
# -------------------------------

cat("Modules created:", length(unique(mod_df$MODULE)), "\n")

table(mod_df$MODULE)
# now run MSEA R script (below) following MSEA Mergeomics tutorial
# https://github.com/XiaYangLabOrg/mergeomics?tab=readme-ov-file#module-merging


# run in bash




# ------------------------------------------------------------------------------

################################################################################
# After wKDA analysis in SCING, create conversion table for cytoscape to convert 
# EnsemblIDs into symbols
################################################################################

# redo for each cell type, changing paths

nodes <- read.delim("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_OPCs/cytoscape/kda2cytoscape.nodes.txt")

genes <- unique(nodes$NODE)
length(genes)

library(biomaRt)

mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")

map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "rgd_symbol"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart
)

library(dplyr)

map <- map %>%
  rename(NodeID = ensembl_gene_id,
         GeneSymbol = external_gene_name)

colnames(map)[1] <- "name"

write.table(
  map,
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_OPCs/cytoscape/rat_symbol_map_network.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# write table for identifying key drivers in cytoscape
# again, redo for all cell types

kda <- read.delim("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_OPCs/cytoscape/kda2cytoscape.top.kds.txt")

key_drivers <- unique(kda$NODNAMES)

write.table(
  key_drivers,
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_OPCs/cytoscape/key_drivers.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------------------------

################################################################################
# after finishing MSEA and compiling results, plot heatmap
################################################################################

# After completing MSEA, run this to compile results and then graph
# change inputs to graph top 50 rat, top 50 mouse, and remaining rat and remaining mouse

library(data.table)

base_dir <- "/u/scratch/v/vturnbil/GWAS_MSEA_Mouse"

gwas_dirs <- list.dirs(base_dir, recursive = FALSE)

read_msea <- function(gwas_dir, region){
  
  msea_dir <- file.path(gwas_dir, region, "msea")
  
  result_file <- list.files(msea_dir,
                            pattern="\\.results\\.txt$",
                            full.names=TRUE)
  
  if(length(result_file)==0) return(NULL)
  
  df <- fread(result_file[1])
  
  # remove control rows
  df <- df[!grepl("^_ctrl", MODULE)]
  
  # keep only what we need
  df <- df[, .(MODULE, FDR)]
  
  df[, GWAS := basename(gwas_dir)]
  df[, REGION := region]
  
  return(df)
}

library(data.table)


# combine MS
ms_results <- rbindlist(lapply(gwas_dirs, read_msea, region="MS"),
                        fill=TRUE)

# save files
fwrite(ms_results, "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/MouseMS_MSEA_FDR.tsv", sep="\t")


library(data.table)
library(ggplot2)
library(tidyr)

ms <- fread("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/MS_MSEA_FDR.tsv")

ms_results <- hp_results[FDR < 0.05]

ms[, logFDR := -log10(FDR)]

ms_matrix <- dcast(ms, GWAS ~ MODULE, value.var="logFDR")

library(tibble)

logFDR_mat <- ms_matrix %>%
  column_to_rownames("GWAS") %>%
  as.matrix()

logFDR_mat <- t(logFDR_mat)

library(ComplexHeatmap)
library(circlize)
library(grid)

ht <- Heatmap(
  logFDR_mat,
  name = "-log10(FDR)",
  
  #################################
  # Color scale
  #################################
  col = colorRamp2(c(0,2,5), c("white","orange","red")),
  
  #################################
  # Explicit heatmap size
  #################################
  width  = unit(ncol(logFDR_mat) * 1.2, "cm"),
  height = unit(nrow(logFDR_mat) * 0.7, "cm"),
  
  #################################
  # Dendrogram
  #################################
  column_dend_height = unit(3, "cm"),
  row_dend_width     = unit(3, "cm"),
  
  #################################
  # Label settings
  #################################
  row_names_gp = gpar(fontsize = 18),
  column_names_gp = gpar(fontsize = 18),
  
  row_names_max_width = max_text_width(
    rownames(logFDR_mat),
    gp = gpar(fontsize = 18)
  ),
  
  column_names_max_height = max_text_height(
    colnames(logFDR_mat),
    gp = gpar(fontsize = 18)
  ),
  
  #################################
  # Titles
  #################################
  column_title = "GWAS Enrichment in Hippocampus DEGs\n(MSEA −log10 FDR)",
  column_title_gp = gpar(fontsize = 36, fontface = "bold"),
  
  #################################
  # Legend
  #################################
  heatmap_legend_param = list(
    title = "-log10(FDR)",
    title_gp = gpar(fontsize = 30, fontface = "bold"),
    labels_gp = gpar(fontsize = 24),
    at = c(0,2,5)
  ),
  
  cluster_rows = TRUE,
  cluster_columns = TRUE
)

pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/MS_GWAS_MSEA_heatmap.pdf",
    width = 120,
    height = 20)

draw(ht)

dev.off()

# ------------------------------------------------------------------------------

# Libraries

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(tibble)

############################################
# Read MSEA results
############################################

ms <- fread("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/MS_MSEA_FDR.tsv")

#mouse
ms <- fread("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/Mouse_MS_MSEA_FDR.tsv")
  
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)

############################################
# 1. Filter significant results
############################################
ms <- ms[FDR < 0.05]

############################################
# 2. Convert to -log10(FDR)
############################################
ms[, logFDR := -log10(FDR)]

############################################
# 3. Build GWAS × MODULE table
############################################
ms_matrix <- dcast(ms, GWAS ~ MODULE, value.var = "logFDR")

############################################
# 4. Convert to matrix (GWAS = rows)
############################################
logFDR_mat <- ms_matrix |>
  column_to_rownames("GWAS") |>
  as.matrix()

############################################
# 5. Replace missing values
############################################
logFDR_mat[is.na(logFDR_mat)] <- 0

############################################
# 6. Remove GWAS with no signal
############################################
logFDR_mat <- logFDR_mat[rowSums(logFDR_mat) > 0, ]

############################################
# 7. Handle infinite values
############################################
if (any(is.infinite(logFDR_mat))) {
  max_val <- max(logFDR_mat[is.finite(logFDR_mat)])
  logFDR_mat[is.infinite(logFDR_mat)] <- max_val
}

colnames(logFDR_mat) <- colnames(ms_matrix)[-1]

colnames(logFDR_mat) <- colnames(logFDR_mat)

colnames(logFDR_mat) <- gsub("_XYT_vs_XXO$", " | Typical", colnames(logFDR_mat))
colnames(logFDR_mat) <- gsub("_XYT_vs_XXT$", " | SCC", colnames(logFDR_mat))
colnames(logFDR_mat) <- gsub("_XXT_vs_XXO$", " | Gonadal", colnames(logFDR_mat))
colnames(logFDR_mat) <- gsub("_XXF_vs_XYM$", " | Typical", colnames(logFDR_mat))
colnames(logFDR_mat) <- gsub("_XX_vs_XY$", " | SCC", colnames(logFDR_mat))
colnames(logFDR_mat) <- gsub("_Ovary_vs_Testis$", " | Gonadal", colnames(logFDR_mat))

# Convert back to FDR
fdr_mat <- 10^(-logFDR_mat)

# Create star matrix
star_mat <- matrix("", nrow = nrow(fdr_mat), ncol = ncol(fdr_mat))
rownames(star_mat) <- rownames(fdr_mat)
colnames(star_mat) <- colnames(fdr_mat)

star_mat[fdr_mat < 0.05]  <- "*"
star_mat[fdr_mat < 0.01]  <- "**"
star_mat[fdr_mat < 0.001] <- "***"

logFDR_mat <- logFDR_mat[
  order(rowMeans(logFDR_mat, na.rm = TRUE), decreasing = TRUE),
]

logFDR_mat <- logFDR_mat[51:128, ] # edit these numbers to include your chosen subset

#check if a term exists in matrix
"ZorinaLichtenwalter_Migraine" %in% rownames(logFDR_mat) # (optional but useful)


col_group <- ifelse(grepl("Typical", colnames(logFDR_mat)), "Typical",
                    ifelse(grepl("Gonadal", colnames(logFDR_mat)), "Gonadal",
                           ifelse(grepl("SCC", colnames(logFDR_mat)), "SCC", "Other")))

col_group <- factor(col_group, levels = c("Typical", "Gonadal", "SCC"))

# Extract components
cell_type <- sub("\\s*\\|.*", "", colnames(logFDR_mat))

# Order by split first, then cell type
ord <- order(col_group, cell_type)

# Apply to both
logFDR_mat <- logFDR_mat[, ord]
col_group <- col_group[ord]

############################################
# 8. Build vertical heatmap
############################################
ht <- Heatmap(
  logFDR_mat,
  name = "-log10(FDR)",
  
  col = colorRamp2(c(0, 2, 5), c("white", "orange", "red")),
  
  width  = unit(ncol(logFDR_mat) * 2, "cm"),
  height = unit(nrow(logFDR_mat) * 1.3, "cm"),
  
  column_dend_height = unit(4, "cm"),
  row_dend_width     = unit(4, "cm"),
  
  row_names_gp = gpar(fontsize = 24),
  column_names_gp = gpar(fontsize = 26),
  column_names_rot = 45,
  
  row_names_max_width = max_text_width(
    rownames(logFDR_mat),
    gp = gpar(fontsize = 24)
  ),
  
  column_names_max_height = unit(6, "cm"),
  
  column_title = "GWAS Enrichment in Rat MS DEGs (Remaining Conditions)\n(MSEA −log10 FDR)",
  column_title_gp = gpar(fontsize = 42, fontface = "bold"),
  
  heatmap_legend_param = list(
    title = "-log10(FDR)",
    title_gp = gpar(fontsize = 36, fontface = "bold"),
    labels_gp = gpar(fontsize = 28),
    at = c(0, 2, 5)
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    fdr <- 10^(-logFDR_mat[i, j])
    stars <- if (fdr < 0.001) {
      "***"
    } else if (fdr < 0.01) {
      "**"
    } else if (fdr < 0.05) {
      "*"
    } else {
      ""
    }
    grid.text(stars, x, y, gp = gpar(fontsize = 14))
  },
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  # ⭐ THIS IS THE KEY ADDITION
  column_split = col_group,
  column_gap = unit(2, "cm")
)

############################################
# 9. Save figure (vertical)
############################################
pdf(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/MS_GWAS_MSEA_heatmap_vertical.pdf",
  width = 40,
  height = 120
)

draw(ht)

dev.off()


################################################################################
# phenotype enrichment analysis (IMPC)
################################################################################


deg_all2 <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg_all2.RDS")

library(dplyr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(purrr)

deg <- deg_all2 %>%
  dplyr::rename(Effect = Comparison) %>%
  dplyr::mutate(
    Effect = recode(
      Effect,
      "XYT_vs_XXO" = "Typical",
      "XXT_vs_XXO" = "Gonadal",
      "XYT_vs_XXT" = "SCC"
    ),
    Group = paste(CellType, Effect, sep = " | ")
  )
unique(deg$CellType)


markers <- deg %>%
  filter(adj.P.Val < 0.05) %>%
  pull(gene) %>%
  unique()

library(dplyr)

# read in data from mouse phenotypes from IMPC database

impc_df <- read_csv("/u/scratch/v/vturnbil/MousePhenotype/Phenotypes.csv")

impc_sets <- impc_df %>%
  group_by(mp_term_name) %>%
  summarise(genes = list(unique(marker_symbol))) %>%
  deframe()

# convert rat genes into mouse orthologs

library(readr)
# Read in a conversion table generated from GeneOrthology
convert <- read_csv("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/HP/MapMyCells/mammalian_orthologs_20231113.csv")

convert_map <- convert %>%
  dplyr::select(Rat_Symbol, Rat_EnsemblID, Mouse_EnsemblID) %>%
  dplyr::filter(!is.na(Rat_EnsemblID), !is.na(Mouse_EnsemblID))

# map SYMBOL → Rat Ensembl
symbol_map <- convert_map %>%
  dplyr::filter(!is.na(Rat_Symbol))

deg2 <- deg %>%
  mutate(gene_original = gene)

symbol_idx <- deg2$gene %in% symbol_map$Rat_Symbol

deg2$gene[symbol_idx] <- symbol_map$Rat_EnsemblID[
  match(deg2$gene[symbol_idx], symbol_map$Rat_Symbol)
]

deg_mouse <- deg2 %>%
  dplyr::inner_join(
    convert_map,
    by = c("gene" = "Rat_EnsemblID")
  ) %>%
  dplyr::mutate(gene = Mouse_EnsemblID)

deg_mouse <- deg_mouse %>%
  dplyr::distinct(CellType, Effect, gene, .keep_all = TRUE)

library(AnnotationDbi)
library(org.Mm.eg.db)

# get all IMPC genes
all_impc_genes <- unique(unlist(impc_sets))

# map SYMBOL → ENSEMBL
impc_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = all_impc_genes,
  keytype = "SYMBOL",
  columns = c("ENSEMBL")
)

impc_sets_ensembl <- impc_df %>%
  inner_join(impc_map, by = c("marker_symbol" = "SYMBOL")) %>%
  group_by(mp_term_name) %>%
  summarise(genes = list(unique(ENSEMBL))) %>%
  deframe()

# significant DEGs
deg_genes <- deg_mouse %>%
  filter(adj.P.Val < 0.05) %>%
  pull(gene) %>%
  unique()

# background universe (IMPORTANT)
background <- deg_mouse$gene %>% unique()


library(purrr)
library(dplyr)

deg_split <- deg_mouse %>%
  group_split(CellType, Effect)

group_names <- deg_mouse %>%
  distinct(CellType, Effect) %>%
  mutate(Group = paste(CellType, Effect, sep = " | ")) %>%
  pull(Group)

names(deg_split) <- group_names

library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)

# map IMPC symbols → mouse Ensembl
impc_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(impc_df$marker_symbol),
  keytype = "SYMBOL",
  columns = c("ENSEMBL")
) %>%
  filter(!is.na(ENSEMBL)) %>%
  distinct()

# build TERM2GENE
term2gene <- impc_df %>%
  inner_join(impc_map, by = c("marker_symbol" = "SYMBOL")) %>%
  dplyr::select(mp_term_name, ENSEMBL) %>%
  distinct()

# run ORA analsis

ora_by_group <- map(deg_split, function(df) {
  
  # get DEGs for this group
  genes <- df %>%
    filter(adj.P.Val < 0.05) %>%
    pull(gene) %>%
    unique()
  
  # skip if too few genes
  if (length(genes) < 10) return(NULL)
  
  enricher(
    gene = genes,
    TERM2GENE = term2gene,
    universe = deg_mouse$gene %>% unique()
  )
})

# create df of results

library(purrr)
library(dplyr)

ora_df <- map2_dfr(
  ora_by_group,
  names(ora_by_group),
  ~{
    if (is.null(.x)) return(NULL)
    
    as.data.frame(.x) %>%
      mutate(Group = .y)
  }
)

# compute -log10(FDR)

ora_df <- ora_df %>%
  mutate(logFDR = -log10(p.adjust))

top_terms <- ora_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(Description) %>%
  summarise(best_p = min(p.adjust)) %>%
  arrange(best_p) %>%
  pull(Description)

library(tidyr)
library(tibble)

logFDR_mat <- ora_df %>%
  filter(Description %in% top_terms) %>%
  dplyr::select(Description, Group, logFDR) %>%
  pivot_wider(
    names_from = Group,
    values_from = logFDR,
    values_fill = 0
  ) %>%
  column_to_rownames("Description") %>%
  as.matrix()

logFDR_mat[logFDR_mat > 10] <- 10

logFDR_mat <- logFDR_mat[rowSums(logFDR_mat) > 0, ]

column_split <- sapply(
  strsplit(colnames(logFDR_mat), " \\| "),
  `[`, 1
)

library(ComplexHeatmap)
library(circlize)
library(grid)

ht <- Heatmap(
  logFDR_mat,
  name = "-log10(FDR)",
  
  #################################
  # Color scale
  #################################
  col = colorRamp2(c(0, 2, 5), c("white", "orange", "red")),
  
  #################################
  # Size
  #################################
  width  = unit(ncol(logFDR_mat) * 2, "cm"),
  height = unit(nrow(logFDR_mat) * 1.3, "cm"),
  
  #################################
  # Clustering
  #################################
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = column_split,
  
  #################################
  # Dendrogram size
  #################################
  column_dend_height = unit(4, "cm"),
  row_dend_width     = unit(4, "cm"),
  
  #################################
  # Labels
  #################################
  row_names_gp = gpar(fontsize = 24),
  column_names_gp = gpar(fontsize = 26),
  column_names_rot = 45,
  
  #################################
  # Label spacing
  #################################
  row_names_max_width = max_text_width(
    rownames(logFDR_mat),
    gp = gpar(fontsize = 24)
  ),
  
  column_names_max_height = unit(6, "cm"),
  
  #################################
  # Title
  #################################
  column_title = "IMPC Enrichment in DEGs\n(ORA −log10 FDR)",
  column_title_gp = gpar(fontsize = 42, fontface = "bold"),
  
  #################################
  # Legend
  #################################
  heatmap_legend_param = list(
    title = "-log10(FDR)",
    title_gp = gpar(fontsize = 36, fontface = "bold"),
    labels_gp = gpar(fontsize = 28),
    at = c(0, 2, 5)
  ),
  
  #################################
  # Significance stars
  #################################
  cell_fun = function(j, i, x, y, width, height, fill) {
    fdr <- 10^(-logFDR_mat[i, j])
    
    stars <- if (fdr < 0.001) {
      "***"
    } else if (fdr < 0.01) {
      "**"
    } else if (fdr < 0.05) {
      "*"
    } else {
      ""
    }
    
    grid.text(stars, x, y, gp = gpar(fontsize = 14))
  }
)

#------------
# 9. Save PDF
#------------
pdf(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/MS_IMPC_ORA_heatmap.pdf",
  width = 40,
  height = 120
)

draw(ht)

dev.off()

# -----------------------------------------------------------------------------

################################################################################
# run MDF in python to start Mergeomics analysis - second time, on mouse data
################################################################################

# processed GWAS files already

# convert rat genes to human orthologs and then make mod and inf files

deg <- load("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/JS_Mouse_DEGs/DEG_df_RNA_SD.rda")

library(tidyverse)
library(readr)

convert <- read_csv(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/HP/MapMyCells/mammalian_orthologs_20231113.csv",
  show_col_types = FALSE
)

colnames(convert)

# Mouse → Human mapping (symbol-based)
map_mouse_symbol <- convert %>%
  filter(
    !is.na(Mouse_Symbol),
    !is.na(Human_Symbol)
  ) %>%
  distinct(Rat_Symbol, .keep_all = TRUE) %>%
  dplyr::select(
    gene = Mouse_Symbol,
    human_gene = Human_Symbol
  )

map_mouse_ens <- convert %>%
  filter(
    !is.na(Mouse_EnsemblID),
    !is.na(Human_Symbol)
  ) %>%
  distinct(Rat_EnsemblID, .keep_all = TRUE) %>%
  dplyr::select(
    gene = Mouse_EnsemblID,
    human_gene = Human_Symbol
  )

deg_human <- DEG_df %>%
  
  # Try symbol-based mapping first
  left_join(map_mouse_symbol, by = "gene") %>%
  
  # Try Ensembl-based mapping for remaining
  left_join(
    map_mouse_ens,
    by = "gene",
    suffix = c("", "_ens")
  ) %>%
  
  # Final mouse gene symbol
  mutate(
    gene_human = coalesce(human_gene, human_gene_ens)
  ) %>%
  
  dplyr::select(-human_gene, -human_gene_ens)

deg_clean <- deg_human %>%
  
  filter(
    adj.P.Val < 0.05,
    !is.na(gene_human),
    !is.na(cell_type)
  ) %>%
  
  group_by(cell_type, gene_human) %>%
  
  slice_max(
    order_by = abs(logFC),
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup()

library(dplyr)

# -------------------------------
# 1. CLEAN INPUT
# -------------------------------

deg_clean <- deg_clean %>%
  filter(!is.na(gene_human), gene_human != "")

# -------------------------------
# 2. SELECT TOP GENES
# per CellType × Comparison
# -------------------------------

deg_clean <- deg_clean %>%
  filter(comparison %in% c("XXF_vs_XYM", "Ovary_vs_Testis", "XX_vs_XY")) %>%
  arrange(adj.P.Val) %>%
  group_by(cell_type, comparison) %>%
  slice_head(n = 500) %>%
  ungroup()

unique(deg_clean$comparison)
# -------------------------------
# 3. SUMMARY CHECK
# -------------------------------

deg_clean %>%
  dplyr::summarise(
    total_rows = n(),
    n_celltypes = n_distinct(cell_type),
    n_comparisons = n_distinct(comparison),
    n_genes = n_distinct(gene_human)
  )

# -------------------------------
# 4. BUILD MODULE FILE
# modules = CellType × Comparison
# -------------------------------

mod_df <- deg_clean %>%
  mutate(module = paste(cell_type, comparison, sep = "_")) %>%
  dplyr::select(module, gene = gene_human) %>%
  distinct()

# -------------------------------
# 5. CLEAN MODULE NAMES
# -------------------------------

mod_df <- mod_df %>%
  mutate(
    module = module %>%
      gsub(" ", "_", .) %>%
      gsub("[^A-Za-z0-9_]", "", .)
  )

colnames(mod_df) <- c("MODULE","GENE")

# -------------------------------
# 6. FILTER SMALL MODULES
# -------------------------------

mod_df <- mod_df %>%
  group_by(MODULE) %>%
  filter(n() >= 10) %>%
  ungroup()

# -------------------------------
# 7. WRITE MODULE FILE
# -------------------------------

write.table(
  mod_df,
  file = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/mouseDE_human.mod",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# -------------------------------
# 8. BUILD INFO FILE
# -------------------------------

inf_df <- mod_df %>%
  distinct(MODULE) %>%
  mutate(
    source = "DEG",
    description = "MS"
  ) %>%
  dplyr::select(MODULE, source, description)

colnames(inf_df) <- c("MODULE","SOURCE","DESCR")

# -------------------------------
# 9. WRITE INFO FILE
# -------------------------------

write.table(
  inf_df,
  file = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/mouseDE_human.inf",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# -------------------------------
# 10. FINAL CHECK
# -------------------------------

cat("Modules created:", length(unique(mod_df$MODULE)), "\n")

table(mod_df$MODULE)
# now run MSEA R script (below) following MSEA Mergeomics tutorial
# https://github.com/XiaYangLabOrg/mergeomics?tab=readme-ov-file#module-merging

# edited run_msea.R to reflect mouse inf and mod files

# run in bash


################################################################################
# prepare data after QC for SATURN analysis
################################################################################

SCT_mouse <- readRDS("/u/scratch/v/vturnbil/SATURN/SCT_mouse_seurat_harmony_sampleid_filtered_annotated.Rds")



# remove groups I don't want from SCT mouse (to match XXO, XXT, XYT groups in rat)
unique(SCT_mouse$Trisomy)
keep_groups <- c("Non-Trisomy")

SCT_mouse <- subset(
  SCT_mouse,
  subset = !(SexChrom == "XY" & Gonad == "Female")
)

SCT_mouse <- subset(
  SCT_mouse,
  subset = Trisomy %in% keep_groups
)

table(SCT_mouse$SexChrom, SCT_mouse$Gonad)


# convert MyClass_Bcas1 into celltype

colnames(MS@meta.data)[
  colnames(MS@meta.data) == "MyClass_Bcas1"
] <- "celltype"

# convert to 5had (mouse and rat both)

## install myfork, which compatible with latest anndata 
# remotes::install_github("zqfang/MuDataSeurat", force = T)
library(MuDataSeurat)

# (optional) step 1: Slim down a Seurat object. So you get raw counts, lognorm counts
# Make a copy of the RNA assay
rna_assay <- mouse[["RNA"]]

# Remove scale.data if it exists
if ("scale.data" %in% names(rna_assay@layers)) {
  rna_assay@layers[["scale.data"]] <- NULL
}

# Put the assay back
mouse[["RNA"]] <- rna_assay

DefaultAssay(mouse) <- "RNA"

# for seurat v5, need to JoinLayer first
DefaultAssay(mouse) = "RNA"
mouse <- JoinLayers(mouse)
# single modality
MuDataSeurat::WriteH5AD(mouse, "/u/scratch/v/vturnbil/SATURN/MS_mouse.h5ad", assay="RNA")



# oops, need to make cell types the same across both species first
library(MuDataSeurat)
library(Seurat)

mouse <- ReadH5AD(
  "/u/scratch/v/vturnbil/SATURN/MS_mouse.h5ad"
)

unique(mouse$celltype)
unique(MS$MyClass_Bcas1)

# rename labels
mouse$celltype_harmonized <- dplyr::recode(
  mouse$celltype,
  
  "Astrocyte" = "Astrocytes",
  
  "OPC" = "OPCs",
  
  "GABA" = "GABAergic Neurons",
  
  "Oligo" = "Oligodendrocytes",
  
  "Glut" = "Glutamatergic Neurons",
  
  "GABA-Chol" = "Cholinergic Neurons",
  
  "Microglia" = "Microglia",
  
  "Vascular" = "Vascular",

  
  .default = mouse$celltype
)

unique(mouse$celltype_harmonized)

MS$celltype_harmonized <- dplyr::recode(
  MS$MyClass_Bcas1,

  "GABAergic Neurons (1)" = "GABAergic Neurons",
  
  "GABAergic Neurons (2)" = "GABAergic Neurons",
  
  "Immature Oligodendrocytes" = "Oligodendrocytes",
  
  "Mature Oligodendrocytes" = "Oligodendrocytes",

  
  .default = MS$MyClass_Bcas1
)

unique(MS$celltype_harmonized)

# convert to 5had (mouse and rat both)

library(Seurat)
library(MuDataSeurat)

# clean metadata
rat@meta.data <- data.frame(
  lapply(rat@meta.data, as.character),
  row.names = rownames(rat@meta.data),
  stringsAsFactors = FALSE
)

Idents(rat) <- rep("cell", ncol(rat))

# extract counts
counts <- GetAssayData(
  rat,
  assay = "RNA",
  layer = "counts"
)

# create old-style assay
rna.assay <- CreateAssayObject(counts = counts)

# fresh object
rat_clean <- CreateSeuratObject(
  counts = rna.assay,
  meta.data = rat@meta.data
)

# remove reductions/graphs
rat_clean@reductions <- list()
rat_clean@graphs <- list()
rat_clean@neighbors <- list()
rat_clean@commands <- list()

# export
MuDataSeurat::WriteH5AD(
  rat_clean,
  "~/downloads/MS_FCG_rat_clean.h5ad",
  assay = "RNA"
)
################################################################################
# find number of degs in mouse data and rat data
################################################################################

load("/u/project/xyang123/jshin/medial_septum_sct/Results/complete_analysis/DEG_processed/DEG_SD_metacell/DEG_df_RNA_SD.rda")

p_val_adj_threshold <- 0.05
log2FC_threshold <- 0.1

library(dplyr)

DF_sig <- DEG_df %>%
  filter(
    adj.P.Val < p_val_adj_threshold,
    abs(logFC) > log2FC_threshold
  )




library(dplyr)
library(tidyr)

# Create presence/absence table
gene_comp <- DF_sig %>%
  filter(comparison %in% c("XXF_vs_XYM",
                           "Ovary_vs_Testis",
                           "XX_vs_XY")) %>%
  distinct(gene, comparison) %>%
  mutate(flag = TRUE) %>%
  pivot_wider(
    names_from = comparison,
    values_from = flag,
    values_fill = FALSE
  )

# Keep only genes present in Typical
gene_comp_typical <- gene_comp %>%
  filter(XXF_vs_XYM)

# Classify Typical genes
gene_comp_typical <- gene_comp_typical %>%
  mutate(
    group = case_when(
      !Ovary_vs_Testis & !XX_vs_XY ~ "Typical_only",
      Ovary_vs_Testis & !XX_vs_XY ~ "Typical+Gonadal",
      !Ovary_vs_Testis & XX_vs_XY ~ "Typical+SCC",
      Ovary_vs_Testis & XX_vs_XY ~ "Typical+Gonadal+SCC"
    )
  )


gene_comp_typical %>%
  count(group, name = "n_DEGs") %>%
  bind_rows(
    tibble(
      group = "Typical_total",
      n_DEGs = nrow(gene_comp_typical)
    ),
    .
  )


library(dplyr)
library(tidyr)

gene_comp <- deg_all2 %>%
  filter(
    adj.P.Val < 0.05,
    Comparison %in% c(
      "XYT_vs_XXO",
      "XXT_vs_XXO",
      "XYT_vs_XXT"
    )
  ) %>%
  distinct(gene, Comparison) %>%
  mutate(flag = TRUE) %>%
  pivot_wider(
    names_from = Comparison,
    values_from = flag,
    values_fill = FALSE
  )

# Keep only genes present in Typical
gene_comp_typical <- gene_comp %>%
  filter(XYT_vs_XXO)

# Classify Typical genes
gene_comp_typical <- gene_comp_typical %>%
  mutate(
    group = case_when(
      !XXT_vs_XXO & !XYT_vs_XXT ~ "Typical_only",
      XXT_vs_XXO & !XYT_vs_XXT ~ "Typical+Gonadal",
      !XXT_vs_XXO & XYT_vs_XXT ~ "Typical+SCC",
      XXT_vs_XXO & XYT_vs_XXT ~ "Typical+Gonadal+SCC"
    )
  )

counts <- gene_comp_typical %>%
  count(group, name = "n_DEGs")

total <- tibble(
  group = "Typical_total",
  n_DEGs = nrow(gene_comp_typical)
)

bind_rows(total, counts)


################################################################################
# after SATURN analysis, construct metacells, perform DEM analysis and construct 
# river plots
################################################################################

#-------------------------------------------------------------------------------
# convert .5had object from SATURN into seurat object
#-------------------------------------------------------------------------------

library(MuDataSeurat)
library(Seurat)


seu <- ReadH5AD(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/Tools/SATURN/SATURN_output/saturn_results/test256_data_MS_mouse_clean_MS_FCG_rat_clean_org_saturn_seed_0.h5ad"
)

seu

Reductions(seu)

Embeddings(seu, "macrogenes")[1:5,1:5]

seu <- RunUMAP(
  seu,
  reduction = "macrogenes",
  dims = 1:50
)

seu

Idents(seu) <- "labels"

umap <- Embeddings(seu, "umap")

p <- DimPlot(seu, reduction = 
              "umap", group.by = "species",label = FALSE)

p

# ------------------------------------------------------------------------------
# Construct Metacells via metacell
# ------------------------------------------------------------------------------

# BiocManager::install("tanaylab/metacell")
# BiocManager::install(c("impute", "preprocessCore", "GO.db"))
# install.packages(c('WGCNA', 'UCell', 'GeneOverlap'))
# install.packages("curl", repos = "https://jeroen.r-universe.dev")
# conda install conda-forge::r-enrichr
# BiocManager::install(c("WGCNA", "UCell", "GenomicRanges", "GeneOverlap"))
# devtools::install_github('smorabit/hdWGCNA', ref='dev',force = TRUE)

library(metacell)


# check the function in package
ls("package:metacell", pattern="tgconfig_mat")


library(Seurat)
library(hdWGCNA)
library(WGCNA)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu) #default 50 PC
seu <- FindNeighbors(seu, reduction = "macrogenes", dims = 1:50)
seu <- FindClusters(seu)
seu <- RunUMAP(seu, dims = 1:50)

Idents(seu) <- "species"

umap <- Embeddings(seu, "umap")

p <- DimPlot(seu, reduction = 
               "umap", group.by = "species",label = FALSE)

p


# put back metadata from SCT_mouse and MS

length(intersect(colnames(seu), colnames(MS)))
length(intersect(colnames(seu), colnames(SCT_mouse)))

for (col in colnames(MS@meta.data)) {
  
  seu[[col]] <- NA
  
  seu@meta.data[common_cells, col] <-
    as.character(MS@meta.data[common_cells, col])
  
}

# generate macrogene embedding to run Metacells on

DefaultAssay(seu)

nrow(seu)

head(rownames(seu))

macro_mat <- t(Embeddings(seu, "macrogenes"))

dim(macro_mat)
# should be 2000 x cells

seu[["SATURN"]] <- CreateAssayObject(
  counts = macro_mat
)

DefaultAssay(seu) <- "SATURN"

nrow(seu[["SATURN"]])

head(rownames(seu[["SATURN"]]))
tail(rownames(seu[["SATURN"]]))

dim(Embeddings(seu, "macrogenes"))

colnames(Embeddings(seu, "macrogenes"))[1:10]


Reductions(seu)

dim(Embeddings(seu, "macrogenes"))

head(colnames(Embeddings(seu, "macrogenes")))

DefaultAssay(seu)

dim(GetAssayData(seu, assay = "SATURN", layer = "counts"))

# add Combined_Groups metadata from MS and SCT_mouse

table(seu$labels2)[1:10]

table(seu$Combined_Sample)[1:10]

common_ms <- intersect(colnames(seu), colnames(MS))

ms_meta <- MS@meta.data[
  common_ms,
  c(
    "sample",
    "group"
  ),
  drop = FALSE
]

seu <- AddMetaData(
  seu,
  metadata = ms_meta
)

common_mouse <- intersect(colnames(seu), colnames(SCT_mouse))

mouse_meta <- SCT_mouse@meta.data[
  common_mouse,
  c(
    "Sample_Name",
    "Genotype",
    "Gonad"
  ),
  drop = FALSE
]

seu <- AddMetaData(
  seu,
  metadata = mouse_meta
)

library(dplyr)

seu$group2 <- case_when(
  seu$Genotype == "XX" & seu$Gonad == "Ovary"  ~ "XXO",
  seu$Genotype == "XX" & seu$Gonad == "Testes" ~ "XXT",
  seu$Genotype == "XY" & seu$Gonad == "Ovary"  ~ "XYO",
  seu$Genotype == "XY" & seu$Gonad == "Testes" ~ "XYT",
  TRUE ~ NA_character_
)

seu$Final_Group <- ifelse(
  is.na(seu$group2),
  seu$group,
  seu$group2
)

table(seu$Final_Group, useNA = "ifany")

seu$Final_Group <- recode(
  seu$Final_Group,
  XXF = "XXO",
  XXM = "XXT",
  XYM = "XYT"
)

seu$Combined_Sample <- ifelse(
  is.na(seu$Sample_Name),
  seu$sample,
  ifelse(
    is.na(seu$sample),
    seu$Sample_Name,
    paste(seu$Sample_Name, seu$sample, sep = "_")
  )
)

table(is.na(seu$Combined_Sample))

table(is.na(seu$Final_Group))

table(is.na(seu$labels2))

table(seu$species, useNA = "ifany")

length(unique(seu$Combined_Sample))

table(seu$group, useNA = "ifany")

table(seu$group2, useNA = "ifany")

table(seu$Final_Group, useNA = "ifany")

head(
  seu@meta.data[
    is.na(seu$Final_Group),
    c(
      "species",
      "group",
      "group2",
      "Genotype",
      "Gonad",
      "sample",
      "Sample_Name"
    )
  ]
)

library(dplyr)

seu$group2 <- case_when(
  # mouse
  seu$Genotype == "XX"  & seu$Gonad == "Ovary"  ~ "XXO",
  seu$Genotype == "XX"  & seu$Gonad == "Testes" ~ "XXT",
  seu$Genotype == "XY"  & seu$Gonad == "Ovary"  ~ "XYO",
  seu$Genotype == "XY"  & seu$Gonad == "Testes" ~ "XYT",
  
  # rat
  seu$Genotype == "XXF" ~ "XXO",
  seu$Genotype == "XXM" ~ "XXT",
  seu$Genotype == "XYM" ~ "XYT",
  
  TRUE ~ NA_character_
)

seu$Final_Group <- dplyr::coalesce(
  seu$group2,
  seu$group
)

table(seu$Final_Group, useNA = "ifany")

table(seu$species, seu$Final_Group)

table(is.na(seu$Combined_Sample))
table(is.na(seu$labels2))
table(is.na(seu$Final_Group))

# get ready for WGCNA

VariableFeatures(seu) <- rownames(seu)

# 1. Setting up the hdWGCNA experiment
# This step is the entry point for hdWGCNA; it adds an empty hdWGCNA experiment to your Seurat object.
seurat_obj <- SetupForWGCNA(
  seu,
  gene_select = "variable", # or "fraction", "topN", "variable"
  #  fraction = 0.05,          # If gene_select = "fraction", select the percentage 
  # of expressed genes. If your goal is to cluster and DEG only on major highly 
  # expressed genes, and the cell number is large, 5% can be used as an initial attempt.
  #To achieve more comprehensive coverage of subgroup markers, it is recommended 
  # to increase the fraction (e.g., 0.10~0.20), or directly use "variable" 
  # screening (e.g., for hypervariable genes).
  # GZ used "variable"
  wgcna_name = "Metacell" # Give your WGCNA experiment a name
)

# varialbe has 2000 genes.
# 0.05 fraction has 4689
length(seurat_obj@misc$Metacell$wgcna_genes)

parallel::detectCores() #check cpu core number

# optionally enable multithreading
enableWGCNAThreads(nThreads = parallel::detectCores())




# 2. Constructing Metacells
# Using MetacellsByGroups()
# If you want to construct Metacells within predefined cell groups (such as Seurat clusters, cell types),
# you can use MetacellsByGroups. This helps reduce sparsity while maintaining cell type specificity.
# For example, constructing metacells within each seurat_clusters:
seurat_obj <- MetacellsByGroups(
  seurat_obj,
  group.by = c("Combined_Sample", "labels2"), # Group by which metadata column 
  # Specify the columns in seurat_obj@meta.data to group by
  k = 10,                       # Number of nearest neighbor cells aggregated within each group's Metacell
  reduction = "macrogenes",
  dims = 1:50,
  layer = "counts",
  mode = "sum",
  min_cells = 10,              # Minimum number of cells required to build a metacell
  max_shared = 1,
  target_metacells = 13272,    # set for ~ 10 cells per metacell
  ident.group ="labels2", # set the Idents of the metacell seurat object
  verbose = TRUE
)



table(seurat_obj$labels2, seurat_obj$Combined_Sample)


# The most convenient way is to directly extract the Seurat object of the 
# Metacell from the hdWGCNA experiment.

metacell_obj <- GetMetacellObject(seurat_obj)


library(dplyr)
library(gridExtra)
library(grid)

# Metacell counts
mc_counts <- as.data.frame(table(metacell_obj$labels2))
colnames(mc_counts) <- c("labels2", "n_metacells")

# Single-cell counts (from Seurat object)
sc_counts <- seu2@meta.data %>%
  dplyr::group_by(labels2) %>%
  dplyr::summarise(n_cells = n(), .groups = "drop")

# Merge + sort
table <- mc_counts %>%
  left_join(sc_counts, by = "labels2") %>%
  arrange(desc(n_metacells))

# Create table grob
table_grob <- tableGrob(
  table,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(cex = 0.8),
      padding = unit(c(2, 2), "mm")
    ),
    colhead = list(fg_params = list(cex = 1))
  )
)


# Save directly to PDF, height scaled to number of rows
pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/celltypes.pdf", width = 6, height = nrow(table)*0.25 + 1)
grid.draw(table_grob)
dev.off()

metacell_obj

# transfer metadata

library(dplyr)

meta_lookup <- seurat_obj@meta.data %>%
  distinct(
    Combined_Sample,
    Final_Group,
    species
  )

meta_lookup %>%
  count(Combined_Sample) %>%
  filter(n > 1)

group_map <- setNames(
  meta_lookup$Final_Group,
  meta_lookup$Combined_Sample
)

species_map <- setNames(
  meta_lookup$species,
  meta_lookup$Combined_Sample
)

metacell_obj@meta.data$Final_Group <-
  unname(group_map[metacell_obj@meta.data$Combined_Sample])

metacell_obj@meta.data$species <-
  unname(species_map[metacell_obj@meta.data$Combined_Sample])

table(is.na(metacell_obj@meta.data$Final_Group))
table(is.na(metacell_obj@meta.data$species))

table(metacell_obj@meta.data$species)

table(is.na(metacell_obj@meta.data$Final_Group))
table(is.na(metacell_obj@meta.data$species))

unique(
  metacell_obj@meta.data$Combined_Sample[
    is.na(metacell_obj@meta.data$Final_Group)
  ]
)

seurat_obj@meta.data %>%
  filter(Combined_Sample %in% unique(
    metacell_obj@meta.data$Combined_Sample[
      is.na(metacell_obj@meta.data$Final_Group)
    ]
  )) %>%
  distinct(Combined_Sample, Final_Group)




seurat_obj@meta.data %>%
  filter(species == "rat") %>%
  select(
    Combined_Sample,
    species,
    group,
    group2,
    Final_Group,
    Genotype,
    Gonad
  ) %>%
  distinct() %>%
  head(30)


seurat_obj$Final_Group <- dplyr::case_when(
  !is.na(seurat_obj$Final_Group) ~ seurat_obj$Final_Group,
  !is.na(seurat_obj$Genotype)    ~ seurat_obj$Genotype
)

table(is.na(seurat_obj$Final_Group))
table(seurat_obj$Final_Group, useNA = "ifany")

seurat_obj$Final_Group <- dplyr::recode(
  seurat_obj$Final_Group,
  "XXF" = "XXO",
  "XXM" = "XXT",
  "XYM" = "XYT"
)

table(seurat_obj$Final_Group, useNA = "ifany")

meta_lookup <- seurat_obj@meta.data %>%
  dplyr::distinct(
    Combined_Sample,
    Final_Group,
    species
  )

group_map <- setNames(
  meta_lookup$Final_Group,
  meta_lookup$Combined_Sample
)

metacell_obj$Final_Group <-
  unname(group_map[metacell_obj$Combined_Sample])

table(metacell_obj$Final_Group, useNA = "ifany")

table(metacell_obj$species, metacell_obj$Final_Group)

colnames(metacell_obj@meta.data)

summary(metacell_obj@meta.data[, c(
  "Combined_Sample",
  "labels2",
  "Final_Group",
  "species"
)])

sapply(
  metacell_obj@meta.data[, c(
    "Combined_Sample",
    "labels2",
    "Final_Group",
    "species"
  )],
  function(x) sum(is.na(x))
)

unique(metacell_obj$Final_Group)


# use the TransferData() from seurat for cell annotation
MSMC <- NormalizeData(metacell_obj)
MSMC <- FindVariableFeatures(MSMC)
MSMC <- ScaleData(MSMC)
MSMC <- RunPCA(MSMC) #default 50 PC
MSMC <- FindNeighbors(MSMC, dims = 1:50)
MSMC <- FindClusters(MSMC)
MSMC <- RunUMAP(MSMC, dims = 1:50)


Idents(MSMC) <- "labels2"

umap <- Embeddings(MSMC, "umap")

p <- DimPlot(MSMC, reduction = 
               "umap", group.by = "labels2",label = TRUE)

p


# Save as MSMC -----------------------------------------------------------------
# saveRDS(MSMC, file = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MSMC_SATURN.RDS")
# MSMC <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Seurat/MSMC_SATURN.RDS")
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#  DEM Analysis with Limma Trend 
# ------------------------------------------------------------------------------

##Script to run limma trend and limma voom for DEG analysis

library(limma)
library(edgeR)
library(data.table)
library(dplyr)
library(plyr)
library(Seurat)
library(org.Mm.eg.db)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggpubr)
library(DESeq2)


##################################################
##load seurat metacell object 
##################################################

library(edgeR)
library(limma)

celltypes <- unique(MSMC$labels2)

deg_list <- list()

assay_to_use <- "SATURN"

for (ct in celltypes) {
  
  for (sp in c("mouse", "rat")) {
    
    message("Processing: ", ct, " | ", sp)
    
    ## ------------------------------------------------------------------
    ## Check whether any metacells exist
    ## ------------------------------------------------------------------
    
    n_cells <- sum(
      MSMC$labels2 == ct &
        MSMC$species == sp &
        MSMC$Final_Group %in% c("XXO", "XXT", "XYT")
    )
    
    if (n_cells == 0) {
      message("  Skipping: no cells found")
      next
    }
    
    ## ------------------------------------------------------------------
    ## Subset
    ## ------------------------------------------------------------------
    
    obj <- subset(
      MSMC,
      subset =
        labels2 == ct &
        species == sp &
        Final_Group %in% c("XXO", "XXT", "XYT")
    )
    
    if (ncol(obj) < 10) {
      message("  Skipping: too few metacells")
      next
    }
    
    meta <- obj@meta.data
    
    meta$Final_Group <- factor(
      meta$Final_Group,
      levels = c("XXO", "XXT", "XYT")
    )
    
    ## ------------------------------------------------------------------
    ## Group counts
    ## ------------------------------------------------------------------
    
    group_counts <- table(meta$Final_Group)
    
    cat("\nGroup counts:\n")
    print(group_counts)
    
    present_groups <- names(group_counts[group_counts >= 3])
    
    if (length(present_groups) < 2) {
      message("  Skipping: fewer than 2 groups with >=3 metacells")
      next
    }
    
    ## Drop unused levels BEFORE model matrix
    meta$Final_Group <- droplevels(meta$Final_Group)
    
    ## ------------------------------------------------------------------
    ## Get counts
    ## ------------------------------------------------------------------
    
    gene_counts <- GetAssayData(
      obj,
      assay = assay_to_use,
      slot = "counts"
    )
    
    gene_counts <- round(gene_counts)
    
    ## ------------------------------------------------------------------
    ## edgeR filtering
    ## ------------------------------------------------------------------
    
    dge0 <- DGEList(gene_counts)
    
    cpm_mat <- edgeR::cpm(dge0)
    
    keep <- apply(cpm_mat, 1, max) >= 2
    
    dge <- dge0[keep, , keep.lib.sizes = FALSE]
    
    cat("\n===== DEG FILTERING SUMMARY =====\n")
    cat("Cell type:", ct, "\n")
    cat("Species:", sp, "\n")
    cat("Metacells:", ncol(obj), "\n")
    cat("Total genes:", nrow(dge0), "\n")
    cat("Genes retained:", nrow(dge), "\n")
    cat("Genes removed:", nrow(dge0) - nrow(dge), "\n")
    cat(
      "Percent retained:",
      round(100 * nrow(dge) / nrow(dge0), 2),
      "%\n"
    )
    cat("=================================\n\n")
    
    if (nrow(dge) < 10) {
      message("  Skipping: too few genes after filtering")
      next
    }
    
    ## ------------------------------------------------------------------
    ## Normalize
    ## ------------------------------------------------------------------
    
    dge <- calcNormFactors(dge)
    
    ## ------------------------------------------------------------------
    ## Design matrix
    ## ------------------------------------------------------------------
    
    design <- model.matrix(
      ~ 0 + Final_Group,
      data = meta
    )
    
    coef_names <- colnames(design)
    
    cat("\nDesign coefficients:\n")
    print(coef_names)
    
    ## Need at least 2 groups
    if (ncol(design) < 2) {
      message("  Skipping: only one group present")
      next
    }
    
    ## ------------------------------------------------------------------
    ## Build valid contrasts dynamically
    ## ------------------------------------------------------------------
    
    contrast_list <- c()
    
    if (all(c("Final_GroupXYT", "Final_GroupXXO") %in% coef_names)) {
      contrast_list <- c(
        contrast_list,
        XYT_vs_XXO = "Final_GroupXYT - Final_GroupXXO"
      )
    }
    
    if (all(c("Final_GroupXYT", "Final_GroupXXT") %in% coef_names)) {
      contrast_list <- c(
        contrast_list,
        XYT_vs_XXT = "Final_GroupXYT - Final_GroupXXT"
      )
    }
    
    if (all(c("Final_GroupXXT", "Final_GroupXXO") %in% coef_names)) {
      contrast_list <- c(
        contrast_list,
        XXT_vs_XXO = "Final_GroupXXT - Final_GroupXXO"
      )
    }
    
    if (length(contrast_list) == 0) {
      message("  Skipping: no valid contrasts")
      next
    }
    
    cat("\nContrasts being tested:\n")
    print(names(contrast_list))
    
    ## ------------------------------------------------------------------
    ## limma
    ## ------------------------------------------------------------------
    
    y <- new("EList")
    
    y$E <- edgeR::cpm(
      dge,
      log = TRUE,
      prior.count = 3
    )
    
    fit <- lmFit(y, design)
    
    contrasts <- do.call(
      makeContrasts,
      c(
        as.list(contrast_list),
        list(levels = design)
      )
    )
    
    fit2 <- contrasts.fit(fit, contrasts)
    
    fit2 <- eBayes(
      fit2,
      trend = TRUE,
      robust = TRUE
    )
    
    ## ------------------------------------------------------------------
    ## Extract results
    ## ------------------------------------------------------------------
    
    resultlist <- list()
    
    for (cmp in colnames(contrasts)) {
      
      tmp <- topTable(
        fit2,
        coef = cmp,
        n = Inf,
        adjust.method = "BH"
      )
      
      tmp$gene <- rownames(tmp)
      tmp$CellType <- ct
      tmp$Species <- sp
      tmp$Comparison <- cmp
      
      resultlist[[cmp]] <- tmp
    }
    
    ## ------------------------------------------------------------------
    ## Store results
    ## ------------------------------------------------------------------
    
    deg_list[[paste(ct, sp, sep = "_")]] <- resultlist
    
    ## ------------------------------------------------------------------
    ## Save per cell type × species
    ## ------------------------------------------------------------------
    
    save(
      resultlist,
      file = paste0(
        "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/",
        gsub(" ", "_", ct),
        "_",
        sp,
        "_pairwise_genotype.rda"
      )
    )
  }
}

## ----------------------------------------------------------------------
## Combine all results
## ----------------------------------------------------------------------

deg_all <- do.call(
  rbind,
  lapply(names(deg_list), function(x) {
    
    do.call(
      rbind,
      deg_list[[x]]
    )
    
  })
)

rownames(deg_all) <- NULL

save(
  deg_all,
  file = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg_all_species.rda"
)

write.csv(
  deg_all,
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg_all_species.csv",
  row.names = FALSE
)


library(dplyr)
library(tidyr)


sig <- deg_all %>%
  dplyr::filter(
    adj.P.Val < 0.05,
    abs(logFC) > 0.1
  ) %>%
  dplyr::mutate(
    direction = ifelse(logFC > 0, "Up", "Down")
  )

mouse_deg <- sig %>%
  dplyr::filter(Species == "mouse") %>%
  dplyr::select(
    CellType,
    Comparison,
    gene,
    mouse_dir = direction
  ) %>%
  dplyr::distinct()

rat_deg <- sig %>%
  dplyr::filter(Species == "rat") %>%
  dplyr::select(
    CellType,
    Comparison,
    gene,
    rat_dir = direction
  ) %>%
  dplyr::distinct()

class(deg_all)
str(deg_all[, 1:10])
colnames(deg_all)


species_class <- dplyr::full_join(
  mouse_deg,
  rat_deg,
  by = c("CellType", "Comparison", "gene")
) %>%
  dplyr::mutate(
    DEG_class = dplyr::case_when(
      !is.na(mouse_dir) & is.na(rat_dir) ~ "Mouse_only",
      
      is.na(mouse_dir) & !is.na(rat_dir) ~ "Rat_only",
      
      !is.na(mouse_dir) &
        !is.na(rat_dir) &
        mouse_dir == rat_dir ~ "Shared_same_direction",
      
      !is.na(mouse_dir) &
        !is.na(rat_dir) &
        mouse_dir != rat_dir ~ "Shared_opposite_direction"
    )
  )

table(species_class$DEG_class)

river_dat <- species_class %>%
  dplyr::group_by(
    CellType,
    Comparison,
    DEG_class
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  )

head(river_dat)

range(river_dat$n)

library(ggplot2)
library(ggalluvial)

ggplot(
  river_dat,
  aes(
    axis1 = CellType,
    axis2 = Comparison,
    y = n
  )
) +
  geom_alluvium(
    aes(fill = DEG_class),
    alpha = 0.8,
    width = 1/12
  ) +
  geom_stratum(
    width = 0.25,
    color = "black"
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      Shared_same_direction     = "#009E73",  # bluish green
      Mouse_only                = "#0072B2",  # blue
      Rat_only                  = "#D55E00",  # vermillion
      Shared_opposite_direction = "#CC79A7"   # reddish purple
    )
  ) +
  scale_x_discrete(
    limits = c(
      "Cell type",
      "Comparison"
    ),
    expand = c(.1, .1)
  ) +
  labs(
    y = "Number of DEMs",
    fill = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12)
  )

ggsave(
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/river_plot.pdf",
  width = 10,
  height = 14
)

# inspect pickle file that has macrogene encodings
# using python script, create csv file with top 50 genes per macrogene

macro_genes <- read.csv("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SATURN/macrogene_top50_genes.csv")

head(macro_genes)

tail(rownames(seu))

deg_all$macrogene <- as.numeric(
  sub("macrogenes-", "", deg_all$gene)
) - 1

sig <- deg_all %>%
  filter(adj.P.Val < 0.05,
         abs(logFC) > 0.1) %>%
  mutate(
    macrogene = as.numeric(
      sub("macrogenes-", "", gene)
    ) - 1
  )

sig_annotated <- sig %>%
  left_join(
    macro_genes,
    by = "macrogene"
  )

sig_annotated %>%
  group_by(
    CellType,
    Comparison,
    macrogene
  ) %>%
  slice_max(weight, n = 10) %>%
  arrange(
    CellType,
    Comparison,
    macrogene,
    desc(weight)
  )

macro_genes <- macro_genes %>%
  tidyr::separate(
    gene,
    into = c("species", "gene_symbol"),
    sep = "_",
    extra = "merge"
  )

head(macro_genes)

table(macro_genes$species)

head(macro_genes[, c("species", "gene_symbol", "macrogene", "weight")])

head(macro_genes$gene)


write.csv(
  sig_annotated,
  file = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SATURN/sig_macrogenes_annotated.csv",
  row.names = FALSE
)







################################################################################
# load JS mouse DEGs
################################################################################

JS <- readRDS("/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/JS_Mouse_DEGs/supercell_seurat_harmony_sampleid_filtered_annotated.Rds")


# remove all but matched genotypes

unique(JS$Genotype)

JS <- subset(JS, Genotype == c("XYM", "XXF", "XXM"))

unique(JS$Genotype)


# supercell has already been run. Do DEG analysis


##Script to run limma trend and limma voom for DEG analysis

library(limma)
library(edgeR)
library(data.table)
library(dplyr)
library(plyr)
library(Seurat)
library(org.Mm.eg.db)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggpubr)
library(DESeq2)


##load seurat metacell object 


celltype_list <- unique(JS$celltype)

deg_list <- list()

assay_to_use <- "RNA"

for (ct in celltype_list) {
  
  message("Processing celltype_list: ", ct)
  
  ## Subset by cell type and genotype
  obj <- subset(
    JS,
    subset = celltype == ct & Genotype %in% c("XXF", "XXM", "XYM")
  )
  
  ## Skip if too few cells
  if (ncol(obj) < 10) {
    message("  Skipping: too few cells")
    next
  }
  
  meta <- obj@meta.data
  meta$Genotype <- factor(meta$Genotype, levels = c("XXF", "XXM", "XYM"))
  meta$Sequencing_Date <- factor(meta$Sequencing_Date)
  
  ## Must have all 3 genotypes for pairwise contrasts
  if (nlevels(droplevels(meta$Genotype)) < 2) {
    message("  Skipping: missing genotypes")
    next
  }
  
  ## Get counts (Seurat v5 compatible)
  gene_counts <- GetAssayData(
    obj,
    assay = assay_to_use,
    slot  = "counts"
  )
  gene_counts <- round(gene_counts)
  
  ## edgeR filtering: keep genes expressed (CPM >= 1) in at least 2 samples
  dge0 <- DGEList(gene_counts)
  cpm_mat <- edgeR::cpm(dge0)
  keep <- apply(cpm_mat, 1, max) >= 2
  dge <- dge0[keep, , keep.lib.sizes = FALSE]
  
  ## ---- Labeled output ----
  cat("\n===== DEG FILTERING SUMMARY =====\n")
  cat("Total genes in assay (before filtering): ",
      nrow(dge0), "\n")
  
  cat("Genes expressed at CPM ≥ 2 in ≥1 sample (tested genes): ",
      nrow(dge), "\n")
  
  cat("Genes removed (low / no expression): ",
      nrow(dge0) - nrow(dge), "\n")
  
  cat("Percent of genes retained for DEG testing: ",
      round(100 * nrow(dge) / nrow(dge0), 2), "%\n")
  
  cat("=================================\n\n")
  
  if (nrow(dge) < 10) {
    message("  Skipping: too few genes after filtering")
    next
  }
  
  dge <- calcNormFactors(dge)
  
  ## Design matrix (try batch correction, remove if confounded)
  
  if (nlevels(meta$Sequencing_Date) > 1) {
    
    design <- model.matrix(
      ~ 0 + Genotype + Sequencing_Date,
      data = meta
    )
    
    if (qr(design)$rank < ncol(design)) {
      
      message(
        "  Sequencing date confounded; removing batch variable"
      )
      
      design <- model.matrix(
        ~ 0 + Genotype,
        data = meta
      )
      
    }
    
  } else {
    
    design <- model.matrix(
      ~ 0 + Genotype,
      data = meta
    )
    
  }
  
  
  ## Final rank check
  if (qr(design)$rank < ncol(design)) {
    
    message(
      "  Skipping: genotype design still rank deficient"
    )
    
    next
    
  }
  
  
  ## Make sure genotype coefficients exist
  required_coef <- c(
    "GenotypeXXF",
    "GenotypeXXM",
    "GenotypeXYM"
  )
  
  if (!all(required_coef %in% colnames(design))) {
    
    message(
      "  Skipping: missing genotype coefficients"
    )
    
    next
    
  }
  
  
  ## Expression matrix for limma
  y <- new("EList")
  y$E <- edgeR::cpm(
    dge,
    log = TRUE,
    prior.count = 3
  )
  
  fit <- lmFit(y, design)
  
  ## Pairwise contrasts
  contrasts <- makeContrasts(
    XYM_vs_XXF = GenotypeXYM - GenotypeXXF,
    XYM_vs_XXM = GenotypeXYM - GenotypeXXM,
    XXM_vs_XXF = GenotypeXXM - GenotypeXXF,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  
  ## Collect results
  resultlist <- list(
    XYM_vs_XXF = topTable(fit2, coef = "XYM_vs_XXF", n = Inf, adjust.method = "BH"),
    XYM_vs_XXM = topTable(fit2, coef = "XYM_vs_XXM", n = Inf, adjust.method = "BH"),
    XXM_vs_XXF = topTable(fit2, coef = "XXM_vs_XXF", n = Inf, adjust.method = "BH")
  )
  
  deg_list[[ct]] <- resultlist
  
  ## Save per cell type
  save(
    resultlist,
    file = paste0(
      "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/",
      gsub(" ", "_", ct),
      "_pairwise_genotype.rda"
    )
  )
}



saveRDS(deg_list, "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/Downsample_Limma_batch.RDS")


deg_list <- readRDS("/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/Downsample_Limma_batch.RDS")


library(dplyr)

# Assume the deg_list structure is deg_list[[ct]][[comp]], and each df row name 
# is either the gene name or the Ensembl ID.
deg_all <- do.call(
  rbind,
  lapply(names(deg_list), function(ct) {
    lapply(names(deg_list[[ct]]), function(comp) {
      df <- deg_list[[ct]][[comp]]
      df$gene <- rownames(df)      # Put the gene name in a separate column
      df$CellType <- ct
      df$Comparison <- comp
      df
    })
  }) %>% unlist(recursive = FALSE)
)

# Row names are currently numbers, and gene names are in the gene column. 
# You can use dplyr::select to adjust the order.
deg_all <- deg_all %>% dplyr::select(gene, CellType, Comparison, everything())

write.csv(deg_all, file = "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/all_DEG_results_Limma.csv", row.names = FALSE)

# plot DEG number

library(readxl)
library(dplyr)
library(ggplot2)

# 2. Set the file path

# Please replace "your_deg_file.xlsx" with your actual Excel file name and path.
file_path <- "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/all_DEG_results_Limma.csv"

# 3. Read Excel file
deg_data <- read.csv(file_path)

# 4. Filter for significant DEGs and count the number of up-regulations/down-regulations

# Define significance threshold and logFC threshold
p_val_adj_threshold <- 0.05    # Adjusted p-value threshold
log2FC_threshold <- 0.1       # avg_log2FC Absolute value greater than 0.1

# Perform filtering and statistics
deg_summary_filtered <- deg_data %>%
  dplyr::filter(adj.P.Val < p_val_adj_threshold) %>%
  dplyr::mutate(
    regulation = dplyr::case_when(
      logFC >  log2FC_threshold  ~ "Up-regulated",
      logFC < -log2FC_threshold  ~ "Down-regulated",
      TRUE                       ~ "Not significant change"
    )
  ) %>%
  dplyr::filter(regulation != "Not significant change") %>%
  dplyr::group_by(CellType, regulation) %>%
  dplyr::summarise(
    count = dplyr::n(),
    .groups = "drop"
  )


# 5. View Results
print(paste0("Statistical results (p_val_adj < ", p_val_adj_threshold, " and |avg_log2FC| > ", log2FC_threshold, "):"))
print(deg_summary_filtered)

# --- Define Color Mapping ---

# Ensure the names here are exactly the same as those you created in the regulation column.

# Here, we'll set Up-regulated to blue and Down-regulated to red, as per your request.

color_mapping <- c("Up-regulated" = "red", "Down-regulated" = "blue")

# If you prefer the standard: Up-regulated red, Down-regulated blue, you can set it like this:

# color_mapping <- c("Up-regulated" = "red", "Down-regulated" = "blue")
# --------------------

# 6. Visualize the results (optional)

# Create a stacked bar chart
p <- ggplot(deg_summary_filtered, aes(x = CellType, y = count, fill = regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = paste0("Number of Up-regulated and Down-regulated DEGs per Cell Type\n(p_val_adj < ", p_val_adj_threshold, " & |avg_log2FC| > ", log2FC_threshold, ")"),
    x = "Cell Type",
    y = "Number of DEGs",
    fill = "Regulation"
  ) +
  scale_y_continuous(expand = c(0, 0)) + # Ensure the y-axis starts from 0 and there is no extra blank space at the top.
  scale_x_discrete(expand = c(0, 0)) +   # Eliminate blank spaces on both sides of the x-axis
  scale_fill_manual(values = color_mapping) + # <-- **Key Modification: Manually Set Colors**
  theme_minimal() + # Basic Theme
  theme(
    panel.grid.major = element_blank(), # Remove main grid lines
    panel.grid.minor = element_blank(), # Remove secondary grid lines
    axis.line = element_line(colour = "black"), # Add coordinate axes and set them to black
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), # Rotate the x-axis labels and set their color.
    axis.text.y = element_text(color = "black"), # Set y-axis label color
    axis.title.x = element_text(color = "black"), # Set x-axis title color
    axis.title.y = element_text(color = "black")  # Set y-axis title color
  )

ggsave(
  filename = "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/Limma_DEGs.pdf",
  plot = p,
  width = 6, height = 4, units = "in",
  dpi = 300
)


################################################################################
# # plot DEG related figures
################################################################################

# Set file path
file_path <- "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/all_DEG_results_Limma.csv"

# 3. Read Excel file
deg_data <- read.csv(file_path)

deg_list <- readRDS("/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/Downsample_Limma_batch.RDS")

library(dplyr)

# Assume the deg_list structure is deg_list[[ct]][[comp]], and each df row name 
# is either the gene name or the Ensembl ID.
deg_all <- do.call(
  rbind,
  lapply(names(deg_list), function(ct) {
    lapply(names(deg_list[[ct]]), function(comp) {
      df <- deg_list[[ct]][[comp]]
      df$gene <- rownames(df)     # Put the gene name in a separate column
      df$CellType <- ct
      df$Comparison <- comp
      df
    })
  }) %>% unlist(recursive = FALSE)
)

# Row names are currently numbers, and gene names are in the gene column. 
# You can use dplyr::select to adjust the order.
deg_all <- deg_all %>% dplyr::select(gene, CellType, Comparison, everything())

rownames(deg_all) <- NULL

# two types of gene names
deg_all <- deg_all %>%
  mutate(gene_type = ifelse(grepl("^ENSMUSG", gene), "Ensembl", "Symbol"))

ens_ids <- unique(deg_all$gene[deg_all$gene_type == "Ensembl"])
symbols <- unique(deg_all$gene[deg_all$gene_type == "Symbol"])

library(biomaRt)
mouse <- useEnsembl(
  biomart = "ensembl",
  dataset = "mmusculus_gene_ensembl",
  mirror = "useast"
)


# Use Ensembl ID to search
ens_annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
  filters = "ensembl_gene_id",
  values = ens_ids,
  mart = mouse
)

# none in ens_ids

# Use symbol to search
symbol_annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
  filters = "external_gene_name",
  values = symbols,
  mart = mouse
)

# combine two tables
gene_annot <- dplyr::bind_rows(
  symbol_annot %>% dplyr::rename(gene = external_gene_name)
)

# no ens_annot

deg_all2 <- deg_all %>%
  left_join(gene_annot, by = "gene") %>%
  mutate(class = ifelse(chromosome_name %in% c("X", "Y"), "Sex-linked", "Autosomal"))

# save deg_all2 ----------------------------------------------------------------
# saveRDS(deg_all2, "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/deg_all2.RDS")
#  deg_all2mouse <- readRDS("/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/deg_all2.RDS")
# ------------------------------------------------------------------------------

write.csv(deg_all2, "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/DEGs.csv")

library(dplyr)
library(ggplot2)
library(forcats)

# 1. Keep only the group XYT_vs_XXO
deg_plot <- deg_all2 %>%
  filter(Comparison == "XYT_vs_XXT" & adj.P.Val < 0.05)

# 2. Group statistics by CellType and class
plotdat <- deg_plot %>%
  dplyr::group_by(CellType, class) %>%
  dplyr::summarise(Freq = n()) %>%
  ungroup()

# 3. Sort CellType in descending order of Freq.
plotdat <- plotdat %>%
  mutate(CellType = fct_reorder(CellType, Freq, .desc = TRUE))
plotdat$facet_label <- "SCC (XYT vs XXT)"

# 4. Draw the drawing and adjust the font size.
p <- ggplot(plotdat, aes(x = CellType, y = log2(Freq), fill = class)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "log2(DEG Count)", fill = "Cell Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), # Adjust font size
        axis.text.y = element_text(angle = 0, hjust = 1, size = 10),
        panel.grid.major = element_blank(),   # Remove main grid lines
        panel.grid.minor = element_blank()  ) + # Remove secondary grid lines
  facet_wrap(~facet_label)

ggsave(
  filename = "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Plots/DEG_plot.pdf",
  plot = p,
  width = 6, height = 4, units = "in",
  dpi = 300
)

#-------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(forcats)

# 1. Keep only the XYT_vs_XXO group, padj < 0.05
deg_plot <- deg_all2 %>%
  filter(Comparison == "XYT_vs_XXT" & adj.P.Val < 0.05)

# 2. Add topping/down tags
deg_plot <- deg_plot %>%
  mutate(Direction = ifelse(logFC > 0, "Up", "Down"))

# 3. Count the number of genes upregulated and downregulated for each cell type
deg_plot <- deg_plot %>%
  dplyr::mutate(Direction = ifelse(logFC > 0, "Up", "Down"))

plotdat <- deg_plot %>%
  dplyr::group_by(CellType, Direction) %>%
  dplyr::summarise(Freq = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(Freq = ifelse(Freq == 0, NA, Freq)) %>%   # avoid log2(0)
  tidyr::pivot_wider(
    names_from = Direction,
    values_from = Freq,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    Up_log2   = ifelse(Up   > 0,  log2(Up),        NA),
    Down_log2 = ifelse(Down > 0, -log2(Down),      NA),
    Total     = Up + Down
  ) %>%
  tidyr::pivot_longer(
    cols = c("Up_log2", "Down_log2"),
    names_to = "Direction",
    values_to = "log2Freq"
  ) %>%
  dplyr::mutate(
    Direction = dplyr::recode(Direction,
                              "Up_log2" = "Up",
                              "Down_log2" = "Down"),
    CellType  = forcats::fct_reorder(CellType, abs(Total), .desc = TRUE),
    facet_label = "SCC (XYT vs XXT)"
  ) %>%
  dplyr::filter(!is.na(log2Freq))

#2. Drawing
p <- ggplot(plotdat, aes(x = CellType, y = log2Freq, fill = Direction)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "log2(DEG Count)", x = NULL, fill = "Direction") +
  facet_wrap(~facet_label) +
  scale_fill_manual(values = c("Down" = "#00bfc4", "Up" = "#f8766d")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/DEG_plot6.pdf",
  plot = p,
  width = 6, height = 4, units = "in",
  dpi = 300
)

# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)

# 1. Keep only the XYM_vs_XXF group
deg_sub <- deg_all2 %>% filter(Comparison == "XX=YM_vs_XXF")

# 2. Extract autosomal and sex-linked data separately.
deg_autosomal <- deg_sub %>% filter(class == "Autosomal")
deg_sexlinked <- deg_sub %>% filter(class == "Sex-linked")

# 3. Aggregate logFC for each (gene, CellType) (take the mean to prevent duplication)
deg_autosomal_unique <- deg_autosomal %>%
  dplyr::group_by(gene, CellType) %>%
  dplyr::summarise(
    logFC = mean(logFC, na.rm = TRUE),
    .groups = "drop"
  )

deg_sexlinked_unique <- deg_sexlinked %>%
  dplyr::group_by(gene, CellType) %>%
  dplyr::summarise(logFC = mean(logFC, na.rm = TRUE), .groups = 'drop')

# 4. Wide table transformation
logFC_mat_autosomal <- deg_autosomal_unique %>%
  pivot_wider(names_from = CellType, values_from = logFC)

logFC_mat_sexlinked <- deg_sexlinked_unique %>%
  pivot_wider(names_from = CellType, values_from = logFC)

# 5. Set row names to numeric matrix
logFC_mat_autosomal <- as.data.frame(logFC_mat_autosomal)
rownames(logFC_mat_autosomal) <- logFC_mat_autosomal$gene
logFC_mat_autosomal$gene <- NULL
logFC_mat_autosomal[] <- lapply(logFC_mat_autosomal, as.numeric)
logFC_mat_autosomal <- as.matrix(logFC_mat_autosomal)

logFC_mat_sexlinked <- as.data.frame(logFC_mat_sexlinked)
rownames(logFC_mat_sexlinked) <- logFC_mat_sexlinked$gene
logFC_mat_sexlinked$gene <- NULL
logFC_mat_sexlinked[] <- lapply(logFC_mat_sexlinked, as.numeric)
logFC_mat_sexlinked <- as.matrix(logFC_mat_sexlinked)

# 6. Correlation Calculation
cor_mat_autosomal <- cor(logFC_mat_autosomal, use = "pairwise.complete.obs", method = "pearson")
cor_mat_sexlinked <- cor(logFC_mat_sexlinked, use = "pairwise.complete.obs", method = "pearson")

cor_mat_autosomal

get_cor_pmat <- function(mat) {
  n <- ncol(mat)
  pmat <- matrix(NA, n, n)
  colnames(pmat) <- rownames(pmat) <- colnames(mat)
  for(i in 1:n) {
    for(j in 1:n) {
      if(i == j) {
        pmat[i, j] <- 1
      } else {
        test <- cor.test(mat[,i], mat[,j], method = "pearson")
        pmat[i, j] <- test$p.value
      }
    }
  }
  pmat
}
p_mat_autosomal <- get_cor_pmat(logFC_mat_autosomal)
p_mat_sexlinked <- get_cor_pmat(logFC_mat_sexlinked)

p_mat_autosomal

library(pheatmap)

# Prepare the asterisk annotation matrix
get_star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}
ann_autosomal <- get_star(p_mat_autosomal)
ann_sexlinked <- get_star(p_mat_sexlinked)

annotation_col = data.frame(
  CellType = colnames(cor_mat_autosomal),
  Effect = effect_vec[colnames(cor_mat_autosomal)]
)
annotation_row = data.frame(
  CellType = rownames(cor_mat_autosomal),
  Effect = effect_vec[rownames(cor_mat_autosomal)]
)

rownames(annotation_col) <- colnames(cor_mat_autosomal)
rownames(annotation_row) <- rownames(cor_mat_autosomal)

unique(annotation_col$CellType)


ann_colors <- list(
  CellType =  c( "Mature Oligodendrocytes" = "darkred",
                 "Immature Oligodendrocytes" = "red",
                 "GABAergic Neurons (1)" = "purple",
                 "GABAergic Neurons (2)" = "lavender",
                 "Astrocytes" = "yellow",
                 "OPCs" = "orange",
                 "Glutamatergic Neurons" = "green",
                 "Cholinergic Neurons" = "blue"
  ), # Your CellType color chart
  Effect = c("Normative" = "purple", "Other" = "lavender") # Custom
)


# heatmap will fail. ignore errors and move on

p <- pheatmap(
  cor_mat_autosomal,
  display_numbers = ann_autosomal,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  main = "Typical (XYT vs XXO) logFC Correlations (Autosomal)",
  fontsize_number = 12,
  breaks = seq(-1, 1, length.out = 101)
)


library(ComplexHeatmap)


pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/Correlation_autosomal.pdf", width=14, height=6) ## ggsave can't save not ggplot
print(p)
dev.off()

# -----------------------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)

# Assuming ann_colors, annotation_col, annotation_row, cor_mat_autosomal,
# and ann_autosomal are already prepared.

ht <- Heatmap(
  cor_mat_autosomal,
  name = "Correlation",  # Main legend title
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  top_annotation = HeatmapAnnotation(
    CellType = annotation_col$CellType,
    col = ann_colors,
    annotation_legend_param = list(
      CellType = list(title = "CellType")
    )
  ),
  left_annotation = rowAnnotation(
    CellType = annotation_row$CellType,
    Effect = annotation_row$Effect,
    col = ann_colors,
    annotation_legend_param = list(
      CellType = list(title = "CellType")
    )
  ),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Gonadal (XXT vs XXO) logFC Correlations
  (Autosomal)",
  heatmap_legend_param = list(
    title = "Correlation",
    at = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    legend_height = unit(2, "cm"),
    border = "black"
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(ann_autosomal[i, j], x, y, gp = gpar(fontsize = 10))
  }
)

# All legends are listed in the column on the right.
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE)

pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/Correlation_autosomal.pdf", width=8, height=7) ## ggsave can't save not ggplot
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

# ------------------------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)
library(grid) # for unit()

# HeatmapAnnotation and rowAnnotation only need to be defined once (consistent with autosomal).）
ht_sexlinked <- Heatmap(
  cor_mat_sexlinked,
  name = "Correlation",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  top_annotation = HeatmapAnnotation(
    CellType = annotation_col$CellType,
    col = ann_colors,
    annotation_legend_param = list(
      CellType = list(title = "CellType")
    )
  ),
  left_annotation = rowAnnotation(
    CellType = annotation_row$CellType,
    Effect = annotation_row$Effect,
    col = ann_colors,
    annotation_legend_param = list(
      CellType = list(title = "CellType")
    )
  ),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Gonadal (XXT vs XXO) logFC Correlations
  (Sex-linked)",
  heatmap_legend_param = list(
    title = "Correlation",
    at = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    legend_height = unit(2, "cm"), # Legend length can be adjusted as needed
    border = "black"
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(ann_sexlinked[i, j], x, y, gp = gpar(fontsize = 10))
  }
)

draw(ht_sexlinked, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE)

pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/Correlation_sex.pdf", width=8, height=7) ## ggsave can't save not ggplot
draw(ht_sexlinked, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()


# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(forcats)

# 1. Only keep 3 comparisons & FDR < 0.05
deg_plot <- deg_all2 %>%
  filter(
    Comparison %in% c("XYT_vs_XXO", "XXT_vs_XXO", "XYT_vs_XXT"),
    adj.P.Val < 0.05
  )

# 2. 设定每个comparison的facet标签，两行显示
comparison_labels <- c(
  "XYT_vs_XXO" = "Typical\n(XYT vs XXO)",
  "XXT_vs_XXO" = "Gonadal\n(XXT vs XXO)",
  "XYT_vs_XXT" = "SCC\n(XYT vs XXT)"
)
deg_plot$facet_label <- comparison_labels[deg_plot$Comparison]

# 3. Group statistics by CellType
plotdat <- deg_plot %>%
  dplyr::group_by(facet_label, CellType) %>%
  dplyr::summarise(Freq = n(), .groups = "drop")

# 4. Set the facet order (Normative, Gonad, Sex chromosome)
plotdat$facet_label <- factor(plotdat$facet_label, levels = c(
  "Typical\n(XYT vs XXO)",
  "Gonadal\n(XXT vs XXO)",
  "SCC\n(XYT vs XXT)"
))

# 5. Set CellType order in descending order of the number of DEG values for Normative features.
normative_order <- plotdat %>%
  filter(facet_label == "Typical\n(XYT vs XXO)") %>%
  arrange(desc(Freq)) %>%
  pull(CellType)
plotdat$CellType <- factor(plotdat$CellType, levels = normative_order)

# 6. Custom Colors
cell_colors <- c(
  "Mature Oligodendrocytes" = "darkred",
  "Immature Oligodendrocytes" = "red",
  "GABAergic Neurons (1)" = "purple",
  "GABAergic Neurons (2)" = "lavender",
  "Astrocytes" = "yellow",
  "OPCs" = "orange",
  "Glutamatergic Neurons" = "green",
  "Cholinergic Neurons" = "blue",
  "Endothelian Cells" = "lightblue",
  "Microglia" = "darkgreen",
  "Vascular" = "gray"
)

# 7. Drawing
p <- ggplot(plotdat, aes(x = CellType, y = log2(Freq), fill = CellType)) +
  geom_bar(stat = "identity", color = "black") +
  theme_bw() +
  labs(y = "log2(DEG Count)", x = "Cell Type", fill = "Cell Type") +
  scale_fill_manual(values = cell_colors, drop = FALSE) +
  scale_y_continuous(expand = c(0, 0), labels = as.integer) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  facet_wrap(~facet_label, nrow = 1, scales = "fixed")

ggsave(
  filename = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_plot8.pdf",
  plot = p,
  width = 8, height = 4, units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(forcats)

# 1. 
deg_plot <- deg_all2 %>%
  filter(
    Comparison %in% c("XYT_vs_XXO", "XXT_vs_XXO", "XYT_vs_XXT"),
    adj.P.Val < 0.05
  ) %>%
  mutate(Direction = ifelse(logFC > 0, "Up", "Down")) 

# 2. 
comparison_labels <- c(
  "XYT_vs_XXO" = "Typical\n(XYT vs XXO)",
  "XXT_vs_XXO" = "Gonadal\n(XXT vs XXO)",
  "XYT_vs_XXT" = "SCC\n(XYT vs XXT)"
)
deg_plot$facet_label <- comparison_labels[deg_plot$Comparison]

# 3. 
plotdat <- deg_plot %>%
  dplyr::group_by(facet_label, CellType, Direction) %>%
  dplyr::summarise(Freq = n(), .groups = "drop")

# 4. 
plotdat$facet_label <- factor(plotdat$facet_label, levels = c(
  "Typical\n(XYT vs XXO)",
  "Gonadal\n(XXT vs XXO)",
  "SCC\n(XYT vs XXT)"
))



# 5. 
normative_order <- plotdat %>%
  dplyr::filter(facet_label == "Typical\n(XYT vs XXO)") %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(
    Total = sum(Freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(Total)) %>%
  dplyr::pull(CellType)

plotdat$CellType <- factor(plotdat$CellType, levels = normative_order)

# 6.
plotdat <- plotdat %>%
  mutate(log2Freq = ifelse(Direction == "Down", -log2(Freq), log2(Freq)))

# 7. 
direction_colors <- c( "Up" = "#f8766d","Down" = "#00bfc4")  # Down=青, Up=粉


# 8. 
ggplot(plotdat, aes(x = CellType, y = log2Freq, fill = Direction)) +
  geom_bar(stat = "identity", color = "black") +
  theme_bw() +
  labs(y = "log2(DEG Count)", x = "Cell Type", fill = "Direction") +
  scale_fill_manual(values = direction_colors) +
  scale_y_continuous(expand = c(0, 0), labels = as.integer) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  facet_wrap(~facet_label, nrow = 1, scales = "fixed")

ggsave("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_all_group_direction.pdf",height = 6,width = 7)

# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(forcats)

# 1. 
deg_plot <- deg_all2 %>%
  filter(
    Comparison %in% c("XYT_vs_XXO", "XXT_vs_XXO", "XYT_vs_XXT"),
    adj.P.Val < 0.05
  )

# 2. 
comparison_labels <- c(
  "XYT_vs_XXO" = "Typical\n(XYT vs XXO)",
  "XXT_vs_XXO" = "Gonadal\n(XXT vs XXO)",
  "XYT_vs_XXT" = "SCC\n(XYT vs XXT)"
)
deg_plot$facet_label <- comparison_labels[deg_plot$Comparison]

# 3. 
plotdat <- deg_plot %>%
  dplyr::group_by(facet_label, CellType, class) %>%
  dplyr::summarise(Freq = n(), .groups = "drop")

# 4.
plotdat$facet_label <- factor(plotdat$facet_label, levels = c(
  "Typical\n(XYT vs XXO)",
  "Gonadal\n(XXT vs XXO)",
  "SCC\n(XYT vs XXT)"
))

# 5. 
normative_order <- plotdat %>%
  dplyr::filter(facet_label == "Typical\n(XYT vs XXO)") %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(Total = sum(Freq)) %>%
  dplyr::arrange(desc(Total)) %>%
  dplyr::pull(CellType)
plotdat$CellType <- factor(plotdat$CellType, levels = normative_order)

# 6. 
plotdat <- plotdat %>%
  mutate(log2Freq = log2(Freq))

# 7.
cell_class_colors <- c("Autosomal" = "#FF6F6F", "Sex-linked" = "#10CFC9")

# 8. 
ggplot(plotdat, aes(x = CellType, y = log2Freq, fill = class)) +
  geom_bar(stat = "identity", color = "black", position = "stack") +
  theme_bw() +
  labs(y = "log2(DEG Count)", x = "Cell Type", fill = "Cell Class") +
  scale_fill_manual(values = cell_class_colors) +
  scale_y_continuous(expand = c(0, 0), labels = as.integer) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  facet_wrap(~facet_label, nrow = 1, scales = "fixed")

ggsave("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_all_group_sexauto.pdf",height = 6,width = 7)

# ------------------------------------------------------------------------------

# # DEG across factors

library(dplyr)
deg_sig <- deg_all2 %>% 
  filter(adj.P.Val < 0.05) %>% 
  dplyr::select(gene, CellType, Comparison,class) # need to add dplyr::, otherwise error

library(dplyr)

dup_counts <- deg_sig %>%
  group_by(gene, CellType, Comparison,class) %>%
  tally() %>%
  filter(n > 1)

# View duplicated combinations and how many times they appear
print(dup_counts)
# cause of your earlier pivot_wider error.tidyr::pivot_wider()`: ! 
# Can't convert `fill` <logical> to <list>. Duplicates must be resolved before
# reshaping the data.

library(dplyr)

deg_sig_nodup <- deg_sig %>%
  distinct(gene, CellType, Comparison,class, .keep_all = TRUE)

sum(duplicated(deg_sig_nodup[, c("gene", "CellType", "Comparison","class")]))
# Should be 0

library(tidyr)
deg_wide <- deg_sig_nodup %>%
  mutate(flag = TRUE) %>%
  pivot_wider(names_from = Comparison, values_from = flag, values_fill = FALSE)

# define genes category
deg_wide <- deg_wide %>%
  mutate(
    deg_class = case_when(
      XYT_vs_XXO & !XXT_vs_XXO & !XYT_vs_XXT ~ "Typical_unique",
      !XYT_vs_XXO & XXT_vs_XXO & !XYT_vs_XXT ~ "Gonadal",
      !XYT_vs_XXO & !XXT_vs_XXO & XYT_vs_XXT ~ "SCC",
      (XYT_vs_XXO + XXT_vs_XXO + XYT_vs_XXT) > 1 ~ "Shared"
    )
  )

library(dplyr)


deg_stats <- deg_wide %>%
  dplyr::group_by(CellType, class, deg_class) %>%
  dplyr::summarise(gene_count = n(), .groups = "drop") %>%
  dplyr::group_by(CellType, class) %>%
  dplyr::mutate(prop = gene_count / sum(gene_count))

deg_stats <- deg_stats %>%
  mutate(deg_class = factor(deg_class, levels = c("Typical_unique", "Shared", "Gonadal", "SCC")))

deg_stats_auto <- deg_stats %>% filter(class == "Autosomal")
deg_stats_sex <- deg_stats %>% filter(class == "Sex-linked")

# Autosomal
order_auto <- deg_stats_auto %>%
  filter(deg_class == "Typical_unique") %>%
  arrange(prop) %>%
  pull(CellType)
# to avoid missing cell type in factors, add cell type without Normative_unique
all_celltypes <- unique(deg_stats_auto$CellType)
order_auto_full <- c(setdiff(all_celltypes, order_auto),order_auto)
deg_stats_auto <- deg_stats_auto %>%
  mutate(CellType = factor(CellType, levels = order_auto_full))

# Sex-linked
order_sex <- deg_stats_sex %>%
  filter(deg_class == "Typical_unique") %>%
  arrange(prop) %>%
  pull(CellType)
# to avoid missing cell type in factors, add cell type without Normative_unique

all_celltypes <- unique(deg_stats_sex$CellType)
order_sex_full <- c(setdiff(all_celltypes, order_sex),order_sex)
deg_stats_sex <- deg_stats_sex %>%
  mutate(CellType = factor(CellType, levels = order_sex_full))

order_auto <- sort(unique(as.character(deg_stats_auto$CellType)))
order_sex  <- sort(unique(as.character(deg_stats_sex$CellType)))

deg_stats_auto <- deg_stats_auto %>%
  mutate(CellType = factor(CellType, levels = order_auto))

deg_stats_sex <- deg_stats_sex %>%
  mutate(CellType = factor(CellType, levels = order_sex))

library(ggplot2)

# Autosomal
p1 <- ggplot(deg_stats_auto, aes(x = CellType, y = prop, fill = deg_class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c( "Shared"="#377eb8", "Typical_unique"="#999999","Gonadal"="#e41a1c", "SCC"="#4daf4a"),
    breaks = c("Typical_unique", "Shared", "Gonadal", "SCC")
  ) +
  scale_y_continuous(expand = c(0, 0)) +   
  theme_bw() +
  labs(title = "Autosomal genes", y = "Proportion", fill = "Category")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 10),
    
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, "cm"),         
    #legend.key.height = unit(1, "cm"),     
    #legend.key.width  = unit(1, "cm")       
  )

# Sex-linked
p2 <- ggplot(deg_stats_sex, aes(x = CellType, y = prop, fill = deg_class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c( "Shared"="#377eb8", "Typical_unique"="#999999","Gonadal"="#e41a1c", "SCC"="#4daf4a"),
    breaks = c("Typical_unique", "Shared", "Gonadal", "SCC")
  ) +
  scale_y_continuous(expand = c(0, 0)) +   
  theme_bw() +
  labs(title = "Sex-linked genes", y = "Proportion", fill = "Category")+ # fill used for changing legend title
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 10),
    
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, "cm") )

ggsave("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_comparison_auto_HP.pdf", plot = p1, height = 4,width = 4)

ggsave("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_comparison_sex_HP.pdf", plot = p2, height = 4,width = 4)


deg_wide %>%
  filter(deg_class %in% c("SCC", "Gonadal")) %>%
  count(deg_class, name = "n_DEGs")
# ------------------------------------------------------------------------------

# add one without separation of autosomal and sex-linked genes

library(dplyr)

deg_stats <- deg_wide %>%
  dplyr::group_by(CellType, deg_class) %>%
  dplyr::summarise(gene_count = n(), .groups = "drop") %>%
  dplyr::group_by(CellType) %>%
  dplyr::mutate(prop = gene_count / sum(gene_count))

deg_stats <- deg_stats %>%
  mutate(deg_class = factor(deg_class, levels = c("Typical_unique", "Shared", "Gonadal", "SCC")))

order_both <- deg_stats %>%
  filter(deg_class == "Typical_unique") %>%
  arrange(prop) %>%
  pull(CellType)
# to avoid missing cell type in factors, add cell type without Normative_unique

all_celltypes <- unique(deg_stats$CellType)
order_sex_full <- c(setdiff(all_celltypes, order_both),order_both)
deg_stats <- deg_stats %>%
  mutate(CellType = factor(CellType, levels = order_sex_full))

order <- sort(unique(as.character(deg_stats$CellType)))

deg_stats <- deg_stats %>%
  mutate(CellType = factor(CellType, levels = order_auto))


p3 <- ggplot(deg_stats, aes(x = CellType, y = prop, fill = deg_class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c( "Shared"="#377eb8", "Typical_unique"="#999999","Gonadal"="#e41a1c", "SCC"="#4daf4a"),
    breaks = c("Typical_unique", "Shared", "Gonadal", "SCC")
  ) +
  scale_y_continuous(expand = c(0, 0)) +   
  theme_bw() +
  labs(title = "All genes", y = "Proportion", fill = "Category")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 10),
    
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, "cm"))

ggsave("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/DEG_comparison_wholeChrom_HP.pdf", plot = p3, height = 4,width = 4)

# p3 looks identical to autosomal genes, but the values are actually slightly different

deg_wide %>%
  filter(deg_class %in% c("SCC", "Gonadal", "Typical_unique", "Shared")) %>%
  count(deg_class, name = "n_DEGs")

# ------------------------------------------------------------------------------
# Identify top 10 up and down-regulated genes per cell type
# ------------------------------------------------------------------------------

deg_all2 <- readRDS("/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Limma/deg_all2.RDS")

deg_sig <- deg_all2 %>%
  dplyr::filter(
    !is.na(logFC),
    !is.na(adj.P.Val),
    adj.P.Val < 0.05,
    Comparison %in% c("XYT_vs_XXO", "XYT_vs_XXT", "XXT_vs_XXO")
  )

top10_up <- deg_sig %>%
  dplyr::filter(logFC > 0) %>%
  dplyr::group_by(CellType, Comparison) %>%
  dplyr::arrange(desc(logFC), .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

top10_down <- deg_sig %>%
  dplyr::filter(logFC < 0) %>%
  dplyr::group_by(CellType, Comparison) %>%
  dplyr::arrange(logFC, .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

top10_degs <- dplyr::bind_rows(
  top10_up   %>% dplyr::mutate(Direction = "Up"),
  top10_down %>% dplyr::mutate(Direction = "Down")
)

table(top10_degs$CellType,top10_degs$Comparison, top10_degs$Direction)

write.csv(
  top10_degs,
  "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/Top10_UpDown_DEGs_by_CellType_by_Comparison.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# make volcanoe plots
# ------------------------------------------------------------------------------


deg_all2 <- readRDS("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg_all2.RDS")


library(ggplot2)
library(dplyr)
library(ggrepel)

volcano_plot <- function(df, celltype, comparison,
                         fdr_cutoff = 0.05,
                         lfc_cutoff = 0.25,
                         top_n = 10) {
  
  df_sub <- df %>%
    filter(
      CellType == celltype,
      Comparison == comparison,
      !is.na(logFC),
      !is.na(adj.P.Val)
    ) %>%
    mutate(
      negLogFDR = -log10(adj.P.Val),
      Significance = case_when(
        adj.P.Val < fdr_cutoff & logFC >= lfc_cutoff  ~ "Up",
        adj.P.Val < fdr_cutoff & logFC <= -lfc_cutoff ~ "Down",
        TRUE                                   ~ "NS"
      )
    )
  
  # select top genes for labeling
  top_genes <- df_sub %>%
    filter(Significance != "NS") %>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n = top_n)
  
  ggplot(df_sub, aes(logFC, negLogFDR)) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c(
      "Up" = "#D62728",
      "Down" = "#1F77B4",
      "NS" = "grey70"
    )) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed") +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf
    ) +
    labs(
      title = paste(celltype, comparison),
      x = "log2 Fold Change",
      y = "-log10(FDR)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

celltypes <- unique(deg_all2$CellType)
comparisons <- c("XYT_vs_XXO", "XYT_vs_XXT", "XXT_vs_XXO")

pdf("/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Plots/VolcanoPlots_All_CellTypes.pdf", width = 6, height = 6)

for (ct in celltypes) {
  for (comp in comparisons) {
    if (nrow(filter(deg_all2, CellType == ct, Comparison == comp)) > 0) {
      print(volcano_plot(deg_all2, ct, comp))
    }
  }
}

dev.off()



################################################################################
# clustering DEGs by correlation analysis
################################################################################

deg_all2 <- readRDS("/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/deg_all2.RDS")

# ==============================================================================
# Libraries
# ==============================================================================
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ==============================================================================
# 1. Rename & relabel comparison → Effect
# ==============================================================================
deg <- deg_all2 %>%
  dplyr::rename(Effect = Comparison) %>%
  dplyr::mutate(
    Effect = recode(
      Effect,
      "XYM_vs_XXF" = "Typical",
      "XXM_vs_XXF" = "Gonadal",
      "XYM_vs_XXM" = "SCC"
    )
  )

# ==============================================================================
# 2. Aggregate logFC per gene × CellType × Effect
# ==============================================================================
deg_unique <- deg %>%
  dplyr::group_by(gene, CellType, Effect) %>%
  dplyr::summarise(
    logFC = mean(logFC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Group = paste(CellType, Effect, sep = " | ")
  )

# ==============================================================================
# 3. Wide matrix: genes × (CellType | Effect)
# ==============================================================================
logFC_mat <- deg_unique %>%
  dplyr::select(gene, Group, logFC) %>%
  pivot_wider(names_from = Group, values_from = logFC)

logFC_mat <- as.data.frame(logFC_mat)
rownames(logFC_mat) <- logFC_mat$gene
logFC_mat$gene <- NULL
logFC_mat[] <- lapply(logFC_mat, as.numeric)
logFC_mat <- as.matrix(logFC_mat)

# ==============================================================================
# 4. Correlation matrix (shared genes only)
# ==============================================================================
cor_mat <- cor(
  logFC_mat,
  use = "pairwise.complete.obs",
  method = "pearson"
)

# ==============================================================================
# 5. Correlation p-value matrix
# ==============================================================================
get_cor_pmat <- function(mat) {
  n <- ncol(mat)
  pmat <- matrix(NA, n, n)
  colnames(pmat) <- rownames(pmat) <- colnames(mat)
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) {
        pmat[i, j] <- 1
      } else {
        pmat[i, j] <- cor.test(mat[, i], mat[, j])$p.value
      }
    }
  }
  pmat
}

p_mat <- get_cor_pmat(logFC_mat)

# ==============================================================================
# 6. Significance stars
# ==============================================================================
get_star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

ann_star <- get_star(p_mat)

# ==============================================================================
# 7. Annotation data (CellType + Effect)
# ==============================================================================
annotation_df <- data.frame(
  Group = colnames(cor_mat)
) %>%
  separate(Group, into = c("CellType", "Effect"), sep = " \\| ")

rownames(annotation_df) <- colnames(cor_mat)
unique(deg_all2$CellType)
# ==============================================================================
# 8. Annotation colors
# ==============================================================================
ann_colors <- list(
  CellType = c(
    "Oligo" = "darkred",
    "Vascular" = "red",
    "GABA" = "purple",
    "GABA-Chol" = "lavender",
    "Astrocyte" = "yellow",
    "OPC" = "orange",
    "Glut" = "green",
    "Microglia" = "blue",
    "Ependymal" = "lightblue"
  ),
  Effect = c(
    "Typical"  = "black",
    "Gonadal"  = "darkgray",
    "SCC"      = "lightgray"
  )
)


# ==============================================================================
# 9. ComplexHeatmap
# ==============================================================================
ht <- Heatmap(
  cor_mat,
  name = "Pearson r",
  
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  
  # -----------------------
  # Dendrogram size (KEY FIX)
  # -----------------------
  column_dend_height = unit(3, "cm"),
  row_dend_width     = unit(3, "cm"),
  
  # Optional: thicker dendrogram lines
  column_dend_gp = gpar(lwd = 3),
  row_dend_gp    = gpar(lwd = 3),
  
  # -----------------------
  # Axis labels
  # -----------------------
  row_names_gp    = gpar(fontsize = 40),
  column_names_gp = gpar(fontsize = 40),
  
  # -----------------------
  # Title
  # -----------------------
  column_title = "logFC Correlations by Cell Type and Effect\n(Autosomal + Sex-linked)",
  column_title_gp = gpar(fontsize = 40, fontface = "bold"),
  
  # -----------------------
  # Heatmap legend
  # -----------------------
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp  = gpar(fontsize = 36, fontface = "bold"),
    labels_gp = gpar(fontsize = 36),
    at = c(-1, 0, 1)
  ),
  
  # -----------------------
  # Top annotation
  # -----------------------
  top_annotation = HeatmapAnnotation(
    df = annotation_df,
    col = ann_colors,
    simple_anno_size = unit(1.2, "cm"),
    annotation_name_gp = gpar(fontsize = 40, fontface = "bold"),
    annotation_legend_param = list(
      title_gp  = gpar(fontsize = 38, fontface = "bold"),
      labels_gp = gpar(fontsize = 36)
    )
  ),
  
  # -----------------------
  # Left annotation
  # -----------------------
  left_annotation = rowAnnotation(
    df = annotation_df,
    col = ann_colors,
    simple_anno_size = unit(1.2, "cm"),
    annotation_name_gp = gpar(fontsize = 40, fontface = "bold"),
    annotation_legend_param = list(
      title_gp  = gpar(fontsize = 38, fontface = "bold"),
      labels_gp = gpar(fontsize = 36)
    )
  ),
  
  show_row_names    = TRUE,
  show_column_names = FALSE,
  
  # -----------------------
  # Cell text (stars)
  # -----------------------
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      ann_star[i, j],
      x, y,
      gp = gpar(fontsize = 20, fontface = "bold")
    )
  }
)


# ==============================================================================
# 10. Save
# ==============================================================================
pdf(
  "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Plots/Correlation_CellType_Effect.pdf",
  width = 42, height = 26
)
draw(
  ht,
  annotation_legend_side = "bottom",
  heatmap_legend_side    = "bottom",
  merge_legend           = TRUE,
  padding = unit(c(20, 5, 10, 250), "mm"))
dev.off()




################################################################################
# orthology analysis mouse vs rat
################################################################################


deg_all2mouse <- readRDS("/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/LimmaJS/deg_all2.RDS")

deg_all2 <- readRDS("/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Limma/deg_all2.RDS")


############################################################
# MOUSE vs RAT DEG OVERLAP
# CELL TYPE × SEX COMPARISON
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)

############################################################
# 1. SETTINGS
############################################################

p_val_adj_threshold <- 0.05
log2FC_threshold <- 0.1

############################################################
# 2. CELL TYPE MAPPING
############################################################

mouse_celltype_map <- c(
  "Astrocyte" = "Astrocytes",
  "Oligo" = "Oligodendrocytes",
  "OPC" = "OPCs",
  "GABA" = "GABAergic neurons",
  "GABA-Chol" = "Cholinergic neurons",
  "Glut" = "Glutamatergic neurons"
)

rat_celltype_map <- c(
  "Astrocytes" = "Astrocytes",
  "Immature Oligodendrocytes" = "Oligodendrocytes",
  "Mature Oligodendrocytes" = "Oligodendrocytes",
  "OPCs" = "OPCs",
  "GABAergic Neurons (1)" = "GABAergic neurons",
  "GABAergic Neurons (2)" = "GABAergic neurons",
  "Cholinergic Neurons" = "Cholinergic neurons",
  "Glutamatergic Neurons" = "Glutamatergic neurons"
)

common_celltypes <- c(
  "Astrocytes",
  "Oligodendrocytes",
  "OPCs",
  "GABAergic neurons",
  "Cholinergic neurons",
  "Glutamatergic neurons"
)

############################################################
# 3. SEX COMPARISON MAPPING
############################################################

# Rat T/O terminology -> Mouse M/F terminology
#
# XYT_vs_XXT -> XYM_vs_XXM
# XYT_vs_XXO -> XYM_vs_XXF
# XXT_vs_XXO -> XXM_vs_XXF

comparison_map <- data.frame(
  Rat_Comparison = c(
    "XYT_vs_XXT",
    "XYT_vs_XXO",
    "XXT_vs_XXO"
  ),
  
  Mouse_Comparison = c(
    "XYM_vs_XXM",
    "XYM_vs_XXF",
    "XXM_vs_XXF"
  ),
  
  Comparison_Label = c(
    "XY vs XX (T)",
    "XY vs XX (O)",
    "XX (T) vs XX (O)"
  ),
  
  stringsAsFactors = FALSE
)

############################################################
# 4. PREPARE MOUSE DATA
############################################################

mouse_df <- deg_all2mouse %>%
  mutate(
    
    # Biological cell type
    Equivalent_CellType =
      unname(mouse_celltype_map[CellType]),
    
    # Keep original mouse comparison
    Mouse_Comparison = Comparison,
    
    # Explicitly identify corresponding rat comparison
    Rat_Comparison = case_when(
      Comparison == "XYM_vs_XXM" ~ "XYT_vs_XXT",
      Comparison == "XYM_vs_XXF" ~ "XYT_vs_XXO",
      Comparison == "XXM_vs_XXF" ~ "XXT_vs_XXO",
      TRUE ~ NA_character_
    ),
    
    # Readable comparison label
    Comparison_Label = case_when(
      Comparison == "XYM_vs_XXM" ~ "XY vs XX (T)",
      Comparison == "XYM_vs_XXF" ~ "XY vs XX (O)",
      Comparison == "XXM_vs_XXF" ~ "XX (T) vs XX (O)",
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    Equivalent_CellType %in% common_celltypes,
    !is.na(Rat_Comparison)
  )

############################################################
# 5. PREPARE RAT DATA
############################################################

rat_df <- deg_all2 %>%
  mutate(
    
    # Biological cell type
    Equivalent_CellType =
      unname(rat_celltype_map[CellType]),
    
    # Keep original rat comparison
    Rat_Comparison = Comparison,
    
    # Explicitly identify corresponding mouse comparison
    Mouse_Comparison = case_when(
      Comparison == "XYT_vs_XXT" ~ "XYM_vs_XXM",
      Comparison == "XYT_vs_XXO" ~ "XYM_vs_XXF",
      Comparison == "XXT_vs_XXO" ~ "XXM_vs_XXF",
      TRUE ~ NA_character_
    ),
    
    # Readable comparison label
    Comparison_Label = case_when(
      Comparison == "XYT_vs_XXT" ~ "XY vs XX (T)",
      Comparison == "XYT_vs_XXO" ~ "XY vs XX (O)",
      Comparison == "XXT_vs_XXO" ~ "XX (T) vs XX (O)",
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    Equivalent_CellType %in% common_celltypes,
    !is.na(Mouse_Comparison)
  )

############################################################
# 6. CHECK COMPARISONS
############################################################

cat("\nMouse comparisons:\n")
print(unique(mouse_df$Mouse_Comparison))

cat("\nRat comparisons:\n")
print(unique(rat_df$Rat_Comparison))

############################################################
# 7. LOAD GENEORTHOLOGY TABLE
############################################################

convert <- read.csv(
  "/u/scratch/v/vturnbil/GSU_FCG/Restart/HP/MapMyCells/mammalian_orthologs_20231113.csv",
  stringsAsFactors = FALSE
)

############################################################
# 8. CLEAN ENSEMBL IDS
############################################################

convert <- convert %>%
  mutate(
    Rat_EnsemblID =
      sub("\\..*$", "", Rat_EnsemblID),
    
    Mouse_EnsemblID =
      sub("\\..*$", "", Mouse_EnsemblID)
  )

############################################################
# 9. CREATE SYMBOL -> ENSEMBL MAPS
############################################################

mouse_symbol_map <- convert %>%
  filter(
    !is.na(Mouse_Symbol),
    !is.na(Mouse_EnsemblID)
  ) %>%
  distinct(
    Mouse_Symbol,
    .keep_all = TRUE
  ) %>%
  select(
    Mouse_Symbol,
    Mouse_EnsemblID
  )

rat_symbol_map <- convert %>%
  filter(
    !is.na(Rat_Symbol),
    !is.na(Rat_EnsemblID)
  ) %>%
  distinct(
    Rat_Symbol,
    .keep_all = TRUE
  ) %>%
  select(
    Rat_Symbol,
    Rat_EnsemblID
  )

############################################################
# 10. MOUSE GENE -> ENSEMBL
############################################################

mouse_to_ensembl <- function(x) {
  
  x_clean <- sub("\\..*$", "", x)
  
  # Already Ensembl
  result <- ifelse(
    grepl("^ENSMUSG", x_clean),
    x_clean,
    NA_character_
  )
  
  # Convert symbols
  idx <- is.na(result)
  
  if (any(idx)) {
    
    result[idx] <-
      mouse_symbol_map$Mouse_EnsemblID[
        match(
          x_clean[idx],
          mouse_symbol_map$Mouse_Symbol
        )
      ]
  }
  
  # Preserve genes that could not be mapped
  idx_unmapped <-
    is.na(result) |
    result == ""
  
  result[idx_unmapped] <-
    paste0(
      "MOUSE_RAW_",
      x_clean[idx_unmapped]
    )
  
  result
}

############################################################
# 11. RAT GENE -> ENSEMBL
############################################################

rat_to_ensembl <- function(x) {
  
  x_clean <- sub("\\..*$", "", x)
  
  # Already Ensembl
  result <- ifelse(
    grepl("^ENSRNOG", x_clean),
    x_clean,
    NA_character_
  )
  
  # Convert symbols
  idx <- is.na(result)
  
  if (any(idx)) {
    
    result[idx] <-
      rat_symbol_map$Rat_EnsemblID[
        match(
          x_clean[idx],
          rat_symbol_map$Rat_Symbol
        )
      ]
  }
  
  # Preserve genes that could not be mapped
  idx_unmapped <-
    is.na(result) |
    result == ""
  
  result[idx_unmapped] <-
    paste0(
      "RAT_RAW_",
      x_clean[idx_unmapped]
    )
  
  result
}

############################################################
# 12. ADD GENE IDS + DEG STATUS
############################################################

mouse_df <- mouse_df %>%
  mutate(
    
    gene_key =
      mouse_to_ensembl(gene),
    
    is_DEG =
      !is.na(adj.P.Val) &
      adj.P.Val < p_val_adj_threshold &
      abs(logFC) > log2FC_threshold
  )

rat_df <- rat_df %>%
  mutate(
    
    gene_key =
      rat_to_ensembl(gene),
    
    is_DEG =
      !is.na(adj.P.Val) &
      adj.P.Val < p_val_adj_threshold &
      abs(logFC) > log2FC_threshold
  )

############################################################
# 13. COLLAPSE DUPLICATE MOUSE GENES
############################################################

mouse_genes <- mouse_df %>%
  group_by(
    Equivalent_CellType,
    Rat_Comparison,
    Mouse_Comparison,
    Comparison_Label,
    gene_key
  ) %>%
  summarise(
    Mouse_DEG = any(is_DEG),
    .groups = "drop"
  )

############################################################
# 14. COLLAPSE DUPLICATE RAT GENES
############################################################

rat_genes <- rat_df %>%
  group_by(
    Equivalent_CellType,
    Rat_Comparison,
    Mouse_Comparison,
    Comparison_Label,
    gene_key
  ) %>%
  summarise(
    Rat_DEG = any(is_DEG),
    .groups = "drop"
  )

############################################################
# 15. TOTAL GENES TESTED
#
# ALL genes tested are retained.
# Orthology is NOT required here.
############################################################

mouse_totals <- mouse_genes %>%
  group_by(
    Equivalent_CellType,
    Rat_Comparison,
    Mouse_Comparison,
    Comparison_Label
  ) %>%
  summarise(
    
    Mouse_total_tested = n(),
    
    Mouse_DEGs =
      sum(Mouse_DEG),
    
    .groups = "drop"
  )

rat_totals <- rat_genes %>%
  group_by(
    Equivalent_CellType,
    Rat_Comparison,
    Mouse_Comparison,
    Comparison_Label
  ) %>%
  summarise(
    
    Rat_total_tested = n(),
    
    Rat_DEGs =
      sum(Rat_DEG),
    
    .groups = "drop"
  )

############################################################
# 16. RAT -> MOUSE ORTHOLOGY
############################################################

rat_mouse_ortholog <- convert %>%
  filter(
    !is.na(Rat_EnsemblID),
    !is.na(Mouse_EnsemblID)
  ) %>%
  select(
    Rat_EnsemblID,
    Mouse_EnsemblID
  ) %>%
  distinct()

############################################################
# 17. CREATE CROSS-SPECIES UNION GENE KEY
############################################################

# Mouse genes use their mouse Ensembl ID.
mouse_genes_union <- mouse_genes %>%
  mutate(
    union_gene = gene_key
  ) %>%
  select(
    Equivalent_CellType,
    Rat_Comparison,
    Mouse_Comparison,
    Comparison_Label,
    union_gene,
    Mouse_DEG
  )

# Rat genes are converted to the corresponding
# mouse Ensembl ID whenever an ortholog exists.
rat_genes_union <- rat_genes %>%
  left_join(
    rat_mouse_ortholog,
    by = c(
      "gene_key" = "Rat_EnsemblID"
    )
  ) %>%
  mutate(
    union_gene = ifelse(
      !is.na(Mouse_EnsemblID),
      Mouse_EnsemblID,
      gene_key
    )
  ) %>%
  select(
    Equivalent_CellType,
    Rat_Comparison,
    Mouse_Comparison,
    Comparison_Label,
    union_gene,
    Rat_DEG
  )

############################################################
# 18. COMBINE MOUSE + RAT
############################################################

combined <- full_join(
  mouse_genes_union,
  rat_genes_union,
  by = c(
    "Equivalent_CellType",
    "Rat_Comparison",
    "Mouse_Comparison",
    "Comparison_Label",
    "union_gene"
  )
) %>%
  mutate(
    Mouse_DEG =
      replace_na(Mouse_DEG, FALSE),
    
    Rat_DEG =
      replace_na(Rat_DEG, FALSE)
  )

############################################################
# 19. FISHER COUNTS
############################################################

fisher_counts <- combined %>%
  group_by(
    Equivalent_CellType,
    Rat_Comparison,
    Mouse_Comparison,
    Comparison_Label
  ) %>%
  summarise(
    
    Shared_DEGs =
      sum(
        Mouse_DEG &
          Rat_DEG
      ),
    
    Mouse_only_DEGs =
      sum(
        Mouse_DEG &
          !Rat_DEG
      ),
    
    Rat_only_DEGs =
      sum(
        !Mouse_DEG &
          Rat_DEG
      ),
    
    Neither_DEG =
      sum(
        !Mouse_DEG &
          !Rat_DEG
      ),
    
    Union_total_tested =
      n(),
    
    .groups = "drop"
  )

############################################################
# 20. ADD TOTALS
############################################################

fisher_results <- fisher_counts %>%
  left_join(
    mouse_totals,
    by = c(
      "Equivalent_CellType",
      "Rat_Comparison",
      "Mouse_Comparison",
      "Comparison_Label"
    )
  ) %>%
  left_join(
    rat_totals,
    by = c(
      "Equivalent_CellType",
      "Rat_Comparison",
      "Mouse_Comparison",
      "Comparison_Label"
    )
  )

############################################################
# 21. FISHER FUNCTION
############################################################

fisher_test_one <- function(
    shared,
    mouse_only,
    rat_only,
    neither
) {
  
  contingency_table <- matrix(
    c(
      shared,
      mouse_only,
      rat_only,
      neither
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  result <- fisher.test(
    contingency_table
  )
  
  data.frame(
    Odds_Ratio =
      unname(result$estimate),
    
    CI_low =
      unname(result$conf.int[1]),
    
    CI_high =
      unname(result$conf.int[2]),
    
    Fisher_P =
      result$p.value
  )
}

############################################################
# 22. RUN FISHER FOR ALL 18 COMBINATIONS
############################################################

fisher_stats <- fisher_results %>%
  rowwise() %>%
  mutate(
    
    Fisher = list(
      fisher_test_one(
        Shared_DEGs,
        Mouse_only_DEGs,
        Rat_only_DEGs,
        Neither_DEG
      )
    )
    
  ) %>%
  unnest(Fisher) %>%
  ungroup()

############################################################
# 23. BH FDR ACROSS ALL 18 TESTS
############################################################

fisher_results <- fisher_stats %>%
  mutate(
    Fisher_FDR =
      p.adjust(
        Fisher_P,
        method = "BH"
      )
  )

############################################################
# 24. CALCULATE OVERLAP
############################################################

fisher_results <- fisher_results %>%
  mutate(
    
    Fold_Overlap = case_when(
      
      Shared_DEGs == 0 ~
        0,
      
      Mouse_DEGs == 0 |
        Rat_DEGs == 0 ~
        NA_real_,
      
      TRUE ~
        (
          Shared_DEGs *
            Union_total_tested
        ) /
        (
          Mouse_DEGs *
            Rat_DEGs
        )
    )
  )

############################################################
# 25. SIGNIFICANCE STARS
############################################################

fisher_results <- fisher_results %>%
  mutate(
    
    significance = case_when(
      
      Fisher_FDR < 0.001 ~ "***",
      
      Fisher_FDR < 0.01 ~ "**",
      
      Fisher_FDR < 0.05 ~ "*",
      
      TRUE ~ ""
    )
  )

############################################################
# 26. ORDER RESULTS
############################################################

fisher_results <- fisher_results %>%
  mutate(
    
    Equivalent_CellType =
      factor(
        Equivalent_CellType,
        levels = common_celltypes
      ),
    
    Comparison_Label =
      factor(
        Comparison_Label,
        levels = c(
          "XY vs XX (T)",
          "XY vs XX (O)",
          "XX (T) vs XX (O)"
        )
      )
  ) %>%
  arrange(
    Equivalent_CellType,
    Comparison_Label
  )

############################################################
# 27. VIEW RESULTS
############################################################

fisher_results %>%
  select(
    Equivalent_CellType,
    Comparison_Label,
    Mouse_Comparison,
    Rat_Comparison,
    Mouse_total_tested,
    Rat_total_tested,
    Union_total_tested,
    Mouse_DEGs,
    Rat_DEGs,
    Shared_DEGs,
    Mouse_only_DEGs,
    Rat_only_DEGs,
    Odds_Ratio,
    Fisher_P,
    Fisher_FDR,
    Fold_Overlap
  )

############################################################
# 28. SAVE RESULTS
############################################################

write.csv(
  fisher_results,
  "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Plots/mouse_rat_Fisher_CellType_Sex_results.csv",
  row.names = FALSE
)

############################################################
# 29. CREATE HEATMAP LABEL
############################################################

heatmap_df <- fisher_results %>%
  mutate(
    Comparison_Name = case_when(
      Rat_Comparison == "XYT_vs_XXO" ~ "Typical",
      Rat_Comparison == "XXT_vs_XXO" ~ "Gonadal",
      Rat_Comparison == "XYT_vs_XXT" ~ "SCC",
      TRUE ~ as.character(Comparison_Label)
    ),
    
    Fold_Overlap_plot = ifelse(
      Shared_DEGs == 0,
      NA_real_,
      Fold_Overlap
    ),
    
    Overlap_Category = ifelse(
      Shared_DEGs == 0,
      "Zero overlap",
      "Overlap"
    ),
    
    FDR_stars = significance
  ) %>%
  mutate(
    Comparison_Name = factor(
      Comparison_Name,
      levels = c(
        "Typical",
        "Gonadal",
        "SCC"
      )
    ),
    
    Equivalent_CellType = factor(
      Equivalent_CellType,
      levels = common_celltypes
    )
  ) %>%
  arrange(
    Equivalent_CellType,
    Comparison_Name
  ) %>%
  mutate(
    x_pos = row_number()
  )

############################################################
# 30. CELL TYPE LABEL PLOT
############################################################

celltype_labels <- data.frame(
  Equivalent_CellType = common_celltypes,
  x = c(2, 5, 8, 11, 14, 17)
)

label_plot <- ggplot(
  celltype_labels,
  aes(
    x = x,
    y = 0.5,
    label = Equivalent_CellType
  )
) +
  geom_text(
    size = 5,
    fontface = "bold"
  ) +
  scale_x_continuous(
    limits = c(0.5, 18.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  theme_void() +
  theme(
    plot.margin = margin(
      5, 10, 0, 10
    )
  )

############################################################
# 31. MAIN HEATMAP
############################################################

heatmap_plot <- ggplot(
  heatmap_df,
  aes(
    x = x_pos,
    y = 1
  )
) +
  
  # Non-zero overlap
  geom_tile(
    data = heatmap_df %>%
      filter(
        Overlap_Category == "Overlap"
      ),
    aes(
      fill = Fold_Overlap_plot
    ),
    color = "black",
    linewidth = 0.5,
    width = 1
  ) +
  
  # Zero overlap
  geom_tile(
    data = heatmap_df %>%
      filter(
        Overlap_Category == "Zero overlap"
      ),
    aes(
      alpha = "Zero overlap"
    ),
    fill = "grey80",
    color = "black",
    linewidth = 0.5,
    width = 1
  ) +
  
  # Thick box around each cell type
  geom_rect(
    data = data.frame(
      xmin = c(
        0.5,
        3.5,
        6.5,
        9.5,
        12.5,
        15.5
      ),
      xmax = c(
        3.5,
        6.5,
        9.5,
        12.5,
        15.5,
        18.5
      ),
      ymin = 0.5,
      ymax = 1.5
    ),
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
    linewidth = 1.5
  ) +
  
  # FDR stars
  geom_text(
    aes(
      label = FDR_stars
    ),
    size = 6,
    fontface = "bold"
  ) +
  
  # Fold overlap
  scale_fill_gradient2(
    name = "Fold Overlap",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 1
  ) +
  
  # Zero overlap legend
  scale_alpha_manual(
    name = NULL,
    values = c(
      "Zero overlap" = 1
    ),
    guide = guide_legend(
      override.aes = list(
        fill = "grey80",
        alpha = 1,
        color = "black"
      )
    )
  ) +
  
  scale_x_continuous(
    breaks = 1:18,
    labels = rep(
      c(
        "Typical",
        "Gonadal",
        "SCC"
      ),
      times = 6
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      0.5,
      1.5
    ),
    breaks = NULL,
    expand = c(
      0,
      0
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 14
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 11
    ),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(
      0,
      10,
      5,
      10
    )
  )

############################################################
# 32. COMBINE LABELS + HEATMAP
############################################################

library(patchwork)

p <- label_plot / heatmap_plot +
  plot_layout(
    heights = c(
      0.35,
      1
    )
  )

############################################################
# 33. DISPLAY
############################################################

print(p)




################################################################################
# BH correction for PharmOmics results
################################################################################


library(readxl)
library(openxlsx)
library(dplyr)

input_file <- "/u/scratch/v/vturnbil/MS_Drug_Results.xlsx"
output_file <- "/u/scratch/v/vturnbil/MS_Drug_Results_BH.xlsx"

# Get sheet names
sheet_names <- excel_sheets(input_file)

# Create output workbook
wb <- createWorkbook()

for (sheet in sheet_names) {
  
  # Read sheet
  df <- read_excel(input_file, sheet = sheet)
  
  # Only perform correction on sheets containing the required columns
  if (all(c("z_scorepvalue", "Drug_network_name") %in% colnames(df))) {
    
    df <- df %>%
      group_by(Drug_network_name) %>%
      mutate(
        z_scorepvalue_BH = p.adjust(
          z_scorepvalue,
          method = "BH"
        )
      ) %>%
      ungroup()
  }
  
  # Write sheet
  addWorksheet(wb, sheet)
  writeData(wb, sheet, df)
}

# Save
saveWorkbook(
  wb,
  output_file,
  overwrite = TRUE
)



################################################################################
# Add MSEA FDR data to GWAS table
################################################################################


library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

############################################################
# 1. READ FILES
############################################################

summary_df <- read_excel("/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Mergeomics/MSEA/Supplementary Table_2.xlsx")

rat <- read.delim(
  "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Mergeomics/MSEA/MS_MSEA_FDR.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

mouse <- read.delim(
  "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Mergeomics/MSEA/Mouse_MS_MSEA_FDR.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


############################################################
# 2. FUNCTION TO FORMAT FDR RESULTS
############################################################

format_fdr <- function(df, species) {
  
  if (species == "Mouse") {
    
    df <- df %>%
      mutate(
        Comparison = case_when(
          str_detect(MODULE, "XXF_vs_XYM$") ~ "Typical",
          str_detect(MODULE, "Ovary_vs_Testis$") ~ "Gonadal",
          str_detect(MODULE, "XX_vs_XY$") ~ "SCC",
          TRUE ~ NA_character_
        ),
        
        CellType = case_when(
          str_detect(MODULE, "^Astrocyte_") ~ "Astrocytes",
          str_detect(MODULE, "^Oligo_") ~ "Oligodendrocytes",
          str_detect(MODULE, "^OPC_") ~ "OPCs",
          str_detect(MODULE, "^GABA_") ~ "GABAergic Neurons",
          str_detect(MODULE, "^Glut_") ~ "Glutamatergic Neurons",
          TRUE ~ NA_character_
        )
      )
    
  } else {
    
    df <- df %>%
      mutate(
        Comparison = case_when(
          str_detect(MODULE, "XYT_vs_XXO$") ~ "Typical",
          str_detect(MODULE, "XXT_vs_XXO$") ~ "Gonadal",
          str_detect(MODULE, "XYT_vs_XXT$") ~ "SCC",
          TRUE ~ NA_character_
        ),
        
        CellType = case_when(
          str_detect(MODULE, "^Astrocytes_") ~ "Astrocytes",
          str_detect(MODULE, "^OPCs_") ~ "OPCs",
          str_detect(MODULE, "^Cholinergic_Neurons_") ~
            "Cholinergic Neurons",
          str_detect(MODULE, "^GABAergic_Neurons_") ~
            "GABAergic Neurons",
          str_detect(MODULE, "^Glutamatergic_Neurons_") ~
            "Glutamatergic Neurons",
          str_detect(MODULE, "^Immature_Oligodendrocytes_") ~
            "Immature Oligodendrocytes",
          str_detect(MODULE, "^Mature_Oligodendrocytes_") ~
            "Mature Oligodendrocytes",
          TRUE ~ NA_character_
        )
      )
  }
  
  df %>%
    filter(
      !is.na(Comparison),
      !is.na(CellType)
    ) %>%
    mutate(
      FDR_text = paste0(
        CellType,
        " (",
        formatC(FDR, format = "f", digits = 3),
        ")"
      )
    ) %>%
    group_by(GWAS, Comparison) %>%
    summarise(
      FDR_text = paste(FDR_text, collapse = "; "),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Comparison,
      values_from = FDR_text,
      names_glue = paste0(species, "_{Comparison}_FDR")
    )
}


############################################################
# 3. FORMAT MOUSE AND RAT RESULTS
############################################################

mouse_fdr <- format_fdr(
  mouse,
  "Mouse"
)

rat_fdr <- format_fdr(
  rat,
  "Rat"
)


############################################################
# 4. COMBINE MOUSE + RAT FDR RESULTS
############################################################

fdr_df <- full_join(
  mouse_fdr,
  rat_fdr,
  by = "GWAS"
)


############################################################
# 5. ADD FDR RESULTS TO SUMMARY TABLE
############################################################

summary_df <- summary_df %>%
  left_join(
    fdr_df,
    by = c("Label" = "GWAS")
  )


############################################################
# 6. DETERMINE SIGNIFICANCE
############################################################

# A species is considered significant if ANY cell type
# has FDR < 0.05 in ANY of the three comparisons.

mouse_sig <- mouse %>%
  mutate(
    Comparison = case_when(
      str_detect(MODULE, "XXF_vs_XYM$") ~ "Typical",
      str_detect(MODULE, "Ovary_vs_Testis$") ~ "Gonadal",
      str_detect(MODULE, "XX_vs_XY$") ~ "SCC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Comparison)) %>%
  group_by(GWAS) %>%
  summarise(
    Mouse_Significant = any(FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

rat_sig <- rat %>%
  mutate(
    Comparison = case_when(
      str_detect(MODULE, "XYT_vs_XXO$") ~ "Typical",
      str_detect(MODULE, "XXT_vs_XXO$") ~ "Gonadal",
      str_detect(MODULE, "XYT_vs_XXT$") ~ "SCC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Comparison)) %>%
  group_by(GWAS) %>%
  summarise(
    Rat_Significant = any(FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  )


############################################################
# 7. ADD SIGNIFICANT Y/N
############################################################

summary_df <- summary_df %>%
  left_join(mouse_sig, by = c("Label" = "GWAS")) %>%
  left_join(rat_sig, by = c("Label" = "GWAS")) %>%
  mutate(
    `Significant in mouse? (Y/N)` = if_else(
      coalesce(Mouse_Significant, FALSE),
      "Y",
      "N"
    ),
    
    `Significant in rat? (Y/N)` = if_else(
      coalesce(Rat_Significant, FALSE),
      "Y",
      "N"
    )
  )


############################################################
# 8. REMOVE TEMPORARY COLUMNS
############################################################

summary_df <- summary_df %>%
  select(
    -Mouse_Significant,
    -Rat_Significant
  )


############################################################
# 9. PUT COLUMNS IN DESIRED ORDER
############################################################

summary_df <- summary_df %>%
  select(
    Label,
    `Reported Trait`,
    
    `Significant in mouse? (Y/N)`,
    Mouse_Typical_FDR,
    Mouse_Gonadal_FDR,
    Mouse_SCC_FDR,
    
    `Significant in rat? (Y/N)`,
    Rat_Typical_FDR,
    Rat_Gonadal_FDR,
    Rat_SCC_FDR,
    
    `Enriched in FCG Mouse? (Y/N)`,
    `Enriched in FCG Rat? (Y/N)`,
    Citation
  )


############################################################
# 10. WRITE EXCEL FILE
############################################################

write.xlsx(
  summary_df,
  "/u/scratch/v/vturnbil/GSU_FCG/Restart/MS/Mergeomics/MSEA/Supplementary_Table_2_with_FDR.xlsx",
  overwrite = TRUE
)
