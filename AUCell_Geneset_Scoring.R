### Score gene signatures / pathways with AUCell

library(data.table)
library(tidyverse)
library(Seurat)
library(AUCell)

## Loading gmt genesets
load_AUCell_genesets <- function(path, ignore_cols = 1){
  x <- scan(path, what="", sep="\n")
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  for(i in 1:ignore_cols){
    y <- lapply(y, `[`, -1) 
  }
  return(y)
}

## AUCell Scoring. Adapted from Gary Bader lab
AUCell_batch <- function(inp_data, genesets, num_batches = 100) {
  ## Scores a data matrix with AUCell in batches. Idea is to limit memory consumption when
  ## scoring with AUCell
  ## INPUTS:
  ##    inp_data = input data, either a dxn matrix of d features, n samples or a Seurat object
  ##                containing such a matrix
  ##    genesets = named list of character vectors, each consisting of a set of gene symbols
  ##    num_batches = number of batches to run AUCell for. More batches = fewer cells (observations)
  ##                  for each batch used for scoring
  ##    slot = slot to use if using a Seurat object
  ##    assay = assay to use if using a Seurat object
  ## RETURNS:
  ##  either an nxp matrix (samples x scores)
  if (is.matrix(inp_data) || is(inp_data, 'dgCMatrix')) {
    num_cells <- ncol(inp_data)
    batch_size <- ceiling(num_cells/num_batches)
    score_mat <- c()
    print('Running AUCell scoring')
    Sys.time()
    for (i in 1:num_batches) {
      print(paste('batch', i, Sys.time()))
      ind1 <- (i-1)*batch_size + 1
      ind2 <- i*batch_size
      if (ind2 > num_cells) {
        ind2 <- num_cells
      }
      gene_rankings <- AUCell::AUCell_buildRankings(inp_data[,ind1:ind2], plotStats = FALSE)
      score_mat_i <- AUCell::AUCell_calcAUC(geneSets = genesets, rankings = gene_rankings)
      score_mat_i <- t(SummarizedExperiment::assay(score_mat_i, 'AUC'))
      score_mat <- rbind(score_mat, score_mat_i)
      gc(full = TRUE, verbose = TRUE)
    }
    print('Finished Scoring')
    print(Sys.time())
    return(score_mat)
  } else if (class(inp_data) == 'Seurat') {
    print('Running AUCell scoring')
    Sys.time()
    score_mat <- AUCell_batch(inp_data = GetAssayData(inp_data, assay = 'RNA'), 
                              genesets = genesets, 
                              num_batches = num_batches)
    colnames(score_mat) <- paste0(colnames(score_mat), '_AUC')
    inp_data <- AddMetaData(inp_data, as.data.frame(score_mat))
    print('Finished Scoring')
    print(Sys.time())
  }
}

# Load genesets
## ignore first two columns because they are name and description. 
geneset = load_AUCell_genesets("example_geneset.gmt", ignore_cols=2)

# Score and get anno. Using raw counts is perfectly ok for this method
  # Runs in batches to save memory
AUCell_scores <- AUCell_batch(GetAssayData(SeuratObject, assay = 'RNA'), genesets = c(geneset), num_batches=10) 
colnames(AUCell_scores) <- colnames(AUCell_scores) %>% stringr::str_replace_all("_AUC.*", "_AUC")

# Add scores to metadata
SeuratObject <- AddMetaData(SeuratObject, as.data.frame(AUCell_scores))
SeuratObject[[]]
