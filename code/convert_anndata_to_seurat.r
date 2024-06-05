#! /home/luanyz0418/miniconda3/bin/Rscript

proj.dir <- "/home/luanyz0418/project/clinical_cases/maojiao"
data.dir <- file.path(proj.dir, 'data')
pub.data.dir <- file.path(proj.dir, 'public_data')

library(sceasy)
h5ad.file = file.path(pub.data.dir, "Full_obj_log_counts_soupx_v2.h5ad")
out.file = file.path(pub.data.dir, "Full_obj_log_counts_soupx_v2.rds")
sceasy::convertFormat(h5ad.file, from="anndata", to="seurat", outFile=out.file)
