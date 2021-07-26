#!/bin/sh

## Run Docker version of pySCENIC on scAML data.
### As per Nature Protocols directions

# Dicklab ##################################################

# GRN from logcounts of single cell AML data
docker run -it --rm \
    -v /Users/andyzeng/Drive/Laptop/Dick_Lab/AMLhierarchies/1_scAML/pySCENIC/data:/scenicdata \
    aertslab/pyscenic:0.10.0 pyscenic grn \
        --num_workers 6 \
        -o /results/scAML.adjacencies.tsv \
        /data/pvg_logcounts.csv \
        /data/lambert2018.txt
        

## CisTarget (mask dropouts recommended)
pyscenic ctx \
results/scAML.adjacencies.tsv \
data/hg19-500bp-upstream-10species.mc9nr.feather \
data/hg19-tss-centered-5kb-10species.mc9nr.feather \
data/hg19-tss-centered-10kb-10species.mc9nr.feather \
--annotations_fname data/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname scenic_scAML_logcounts.csv \
--mode "dask_multiprocessing" \
--output results/scAML_maskdropout.regulons.csv \
--num_workers 16 \
--mask_dropouts


## AUCell (just on logcounts)
pyscenic aucell \
scenic_scAML_logcounts.csv \
results/scAML_maskdropout.regulons.csv \
--output results/scAML_maskdropout.auc_mtx.csv \
--num_workers 16

