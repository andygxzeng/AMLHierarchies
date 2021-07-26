# # Scran Sam PVG AML
# **R Script to scran normalize PVG AML samples.**
#
# Andy Zeng 

library(tidyverse)
library(data.table)
library(SingleCellExperiment)
library(scater)
library(scran)

library(readr)


# function to load SCE
load_pvg_AML <- function(exprs_file, anno_file){
  exprs <- fread(exprs_file, sep = "\t")
  anno <- read_delim(anno_file, delim = "\t", col_types = cols()) %>% dplyr::rename("index" = Cell)
  
  # remove header row, get counts and genes
  counts <- exprs[which(exprs[,1]!="Gene"),-1]
  genes <- exprs[which(exprs[,1]!="Gene"),1]
  
  # matrix of counts
  count_matrix <- data.matrix(counts) %>% as(., "dgCMatrix")
  rownames(count_matrix) <- genes[[1]]
  sce <- SingleCellExperiment(list(counts = count_matrix))

  # add coldata
  for(i in 1:ncol(anno)){ colData(sce)[paste(colnames(anno)[i])] <- anno[i] }
  colnames(sce) <- anno$index
  colData(sce)$index <- NULL
  
  return(sce)
}


pvg_AML_preprocessing <- function(sce, min_cells = 0, spike_in = FALSE){
  # QC Metrics
  sce <- calculateQCMetrics(sce)

  # Quick Filter: keep genes that are expressed in at least n cells
  num.cells <- nexprs(sce, byrow=TRUE)
  keep.genes <- which(num.cells >= min_cells)
  sce <- sce[keep.genes,]
  
  # Normalize using pooled size factors + deconvolution
  sce <- computeSumFactors(sce)
  if(spike_in == TRUE){
      sce <- computeSpikeFactors(sce)
  }
  sce <- normalize(sce, log_exprs_offset = 1)
  
  return(sce)
}

### Other samples #################################
do_scran_normalization <- function(sampname, dir, dem, anno, min_cells = 5){
  dat = load_pvg_AML(paste0(dir, dem), paste0(dir, anno)) %>%
    pvg_AML_preprocessing(., min_cells = min_cells)
  
  # Print out new GEM and Anno
  logcounts(dat) %>% as.matrix() %>% as.data.frame() %>% 
    rownames_to_column(var = "gene") %>% write_csv(paste0(sampname,"_AML_hierarchy_normalized_gem.csv"))
  colData(dat) %>% as.data.frame() %>% 
    rownames_to_column(var = "Cell") %>% write_csv(paste0(sampname, "_AML_hierarchy_normalized_anno.csv"))
  
  print(dat)
}

# Normalize one at a time for initial cluster assignment
dir = "../../../../../../pvg_AML_samples/data/other/"

# Samples with most HSC-like and Prog-like cells
do_scran_normalization("AML328", dir, "GSM3587931_AML328-D0.dem.txt", "GSM3587932_AML328-D0.anno.txt", 5)
do_scran_normalization("AML916", dir, "GSM3587988_AML916-D0.dem.txt", "GSM3587989_AML916-D0.anno.txt", 5)
do_scran_normalization("AML921A", dir, "GSM3587990_AML921A-D0.dem.txt", "GSM3587991_AML921A-D0.anno.txt", 5)
do_scran_normalization("AML707B_D0", dir, "GSM3587969_AML707B-D0.dem.txt", "GSM3587970_AML707B-D0.anno.txt", 5)

# Other samples
do_scran_normalization("AML1012", dir, "GSM3587923_AML1012-D0.dem.txt", "GSM3587924_AML1012-D0.anno.txt", 5)
do_scran_normalization("AML210A", dir, "GSM3587925_AML210A-D0.dem.txt", "GSM3587926_AML210A-D0.anno.txt", 5)
do_scran_normalization("AML329", dir, "GSM3587940_AML329-D0.dem.txt", "GSM3587941_AML329-D0.anno.txt", 5)
do_scran_normalization("AML419A", dir, "GSM3587950_AML419A-D0.dem.txt", "GSM3587951_AML419A-D0.anno.txt", 5)
do_scran_normalization("AML420B", dir, "GSM3587953_AML420B-D0.dem.txt", "GSM3587954_AML420B-D0.anno.txt", 5)
do_scran_normalization("AML475", dir, "GSM3587959_AML475-D0.dem.txt", "GSM3587960_AML475-D0.anno.txt", 5)
do_scran_normalization("AML870", dir, "GSM3587984_AML870-D0.dem.txt", "GSM3587985_AML870-D0.anno.txt", 5)
do_scran_normalization("AML556", dir, "GSM3587963_AML556-D0.dem.txt", "GSM3587964_AML556-D0.anno.txt", 5)

