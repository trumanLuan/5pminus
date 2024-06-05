
## conda env: r4.1

proj.dir <- "/home/luanyz0418/project/clinical_cases/maojiao"
data.dir <- file.path(proj.dir, 'data')
pub.data.dir <- file.path(proj.dir, 'public_data')

## required packages
require(readr)
require(dplyr)
require(stringr)
require(Seurat)
require(reshape2)
# library(SummarizedExperiment) 

# library(factoextra)
# library("FactoMineR")

library(ComplexHeatmap)



## ----------------------------
## load genome vars in case2 and case1
## --------------------------

snv.case1 <- read.table(file.path(data.dir, 'case1/case1.snv.final.tsv'), header=T, sep='\t', as.is=T )
cnv.case1 <- read.table(file.path(data.dir, 'case1/case1.cnv.final.tsv'), header=T, sep='\t', as.is=T )

snv.case2 <- readr::read_rds( file.path(data.dir, "case2/case2.snvs.filter.rds") ) %>% as.data.frame()
cnv.case2 <- readr::read_rds( file.path(data.dir, "case2/case2.cnvs.filter.rds") )
sv.case2 <- readr::read_rds( file.path(data.dir, "case2/case2.svs.filter.rds") )

## add read depth column into snv.case2 and sv.case2.

tmp.depth <- sapply(snv.case2$INFO, function(x) stringr::str_extract(x, "DP=\\d+"))
tmp.depth <- stringr::str_replace_all(tmp.depth, "DP=", "") %>% as.numeric()
snv.case2$DEPTH = tmp.depth

sv.case2$SU <- sapply(sv.case2$info, function(x) stringr::str_split(x, ":")[[1]][2])
sv.case2$SU <- as.numeric(sv.case2$SU)

snv.case2.new <- snv.case2 %>% filter(DEPTH > 20)
sv.case2.new <- sv.case2 %>% filter(SU > 20)

## exonic and splicing variants related genes.

snv.case1.select <- snv.case1
cnv.case1.select <- cnv.case1 %>% filter(region.final %in% c("exonic", 'intronic', 'splicing', 'UTR3', 'UTR5'))

snv.case2.select <- snv.case2.new %>% filter(`Func.refGene` %in% c("exonic", 'splicing', 'UTR3', 'UTR5'))
cnv.case2.select <- cnv.case2 %>% filter(gene_region %in% c( "exonic", 'intronic', 'splicing', 'UTR3', 'UTR5' ))
sv.case2.select <- sv.case2.new %>% filter(gene_region %in% c( "exonic", 'intronic', 'splicing', 'UTR3', 'UTR5', "UTR5;UTR3" ))

cnv.case2.select.del <- cnv.case2.select %>% filter(svtype == 'loss')
sv.case2.select.del <- sv.case2.select %>% filter(svtype == 'DEL')

gene.case1 <- c(snv.case1.select$Gene.refGene, cnv.case1.select$genes.final)
gene.case2 <- c(snv.case2.select$Gene.refGene, cnv.case2.select.del$genes, sv.case2.select.del$genes)

gene.case1 <- unique(unlist(sapply(gene.case1, function(x) stringr::str_split(x, ";")[[1]]) ))
gene.case2 <- unique(unlist(sapply(gene.case2, function(x) stringr::str_split(x, ';')[[1]]) ))


gene.shared <- Reduce(intersect, list(gene.case1, gene.case2))


## ----------------------------
## load gutcellatlas single cell gene expression data
## --------------------------
## .h5ad file was converted to rds file using script convert_anndata_to_seurat.r

## .....
## load gutcellatlas data
## first trim, 
## second trim, 
## MLN, mesenteric lymph nodes
## 

gutcell.obj <- readr::read_rds(file.path(pub.data.dir, "Full_obj_log_counts_soupx_v2.rds"))

gutcell.data <- gutcell.obj@assays$RNA@data
gutcell.meta <- gutcell.obj@meta.data

## whether gene.case1 and gene.case2 were identifed in gutcell.data
table(gene.case1 %in% rownames( gutcell.data) )
table(gene.case2 %in% rownames( gutcell.data ) )
table(gene.shared %in% rownames( gutcell.data ) ) # 7 out of 8

# library(SeuratDisk)
# loom.file <- file.path(pub.data.dir, "allgutcombined.loom")
# test <- Connect(filename = loom.file, mode='r')
# test$close_all()

# ## error occured.
# h5.file <- file.path(pub.data.dir, "Full_obj_log_counts_soupx_v2.h5ad")
# Convert(h5.file, ".h5seurat") ## error occured. 
# seuratObj <- LoadH5Seurat(file.path(pub.data.dir, "Full_obj_log_counts_soupx_v2.h5Seurat") )


# ## error occured. 
# library(anndata)
# library(reticulate)
# sce = anndata::read_h5ad(file.path(pub.data.dir, "Full_obj_log_counts_soupx_v2.h5ad")) # answer 'no'
# sce.seurat = Convert(sce, to = "seurat")


# ## --------------
# ## load ENCC single cell rnaseq data
# ## --------------

# hscr.data = file.path(pub.data.dir, "HSCR.seurat.Rdata")
# load(hscr.data) #HSCR.seurat
# rm(hscr.data)

# encc.obj <- HSCR.seurat

# encc.meta <- encc.obj@meta.data
# encc.meta$cell <- rownames(encc.meta)

# encc.data <- encc.obj@assays$RNA@data

# ## check overlap between gene vars and genes included in single cell data
# table(gene.case1 %in% rownames( encc.data) )
# table(gene.case2 %in% rownames( encc.data ) )
# table(gene.shared %in% rownames( encc.data ) ) # 7 out of 8


# ## check missing rate of gene expr values
# encc.subset.case1 <- encc.data[ rownames( encc.data ) %in% gene.case1,]
# encc.subset.case2 <- encc.data[ rownames( encc.data ) %in% gene.case2,]
# encc.subset.shared <- encc.data[ rownames( encc.data ) %in% gene.shared,]

# encc.subset.shared <- as.matrix(encc.subset.shared)
# encc.subset.case2 <- as.matrix(encc.subset.case2)
# encc.subset.case1 <- as.matrix(encc.subset.case1)


# # encc.case1.melt <- reshape2::melt(encc.subset.case1)
# # encc.case2.melt <- reshape2::melt(encc.subset.case2)
# # encc.shared.melt <- reshape2::melt(encc.subset.shared)

# # colnames(encc.case1.melt) <- c('gene', 'cell', 'exp')
# # colnames(encc.case2.melt) <- c('gene', 'cell', 'exp')
# # colnames(encc.shared.melt) <- c('gene', 'cell', 'exp')

# # encc.case1.melt$iszero <- encc.case1.melt$exp == 0
# # encc.case2.melt$iszero <- encc.case2.melt$exp == 0
# # encc.shared.melt$iszero <- encc.shared.melt$exp == 0

# # map_cell_id <- function(assay.melt, meta.data){
# # 	t = sapply(assay.melt$cell, function(x) subset(meta.data, cell == x)$severity)
# # 	assay.melt$severity = t
# # 	return(assay.melt)
# # }

# # encc.case1.melt <- map_cell_id(encc.case1.melt, encc.meta)
# # readr::write_rds(encc.case1.melt, file.path(pub.data.dir, "encc.subset.with.case1.genes.rds"))

# # encc.case2.melt <- map_cell_id(encc.case2.melt, encc.meta)
# # readr::write_rds(encc.case2.melt, file.path(pub.data.dir, "encc.subset.with.case2.genes.rds"))

# # encc.shared.melt <- map_cell_id(encc.shared.melt, encc.meta)
# # readr::write_rds(encc.shared.melt, file.path(pub.data.dir, "encc.subset.with.shared.genes.in.two.cases.rds"))

# ## //
# encc.case1.melt <- readr::read_rds( file.path(pub.data.dir, "encc.subset.with.case1.genes.rds") )
# encc.case2.melt <- readr::read_rds( file.path(pub.data.dir, "encc.subset.with.case2.genes.rds") )
# encc.shared.melt <- readr::read_rds( file.path(pub.data.dir, "encc.subset.with.shared.genes.in.two.cases.rds") )

# a <- encc.shared.melt %>% group_by(gene, severity) %>% summarize(zero_count=sum(iszero), total_count=n_distinct(cell)) %>%
# 	as.data.frame()
# a$missing.rate <- a$zero_count / a$total_count

# a.new <- a %>% filter(missing.rate <0.5)


# ## //

# a2 <- encc.case1.melt %>% group_by(gene, severity) %>% summarize(zero_count=sum(iszero), total_count=n_distinct(cell)) %>%
# 	as.data.frame()
# a2$missing.rate <- a2$zero_count / a2$total_count

# a2.new <- a2 %>% filter(missing.rate <0.5)
# a2.new.in <- a2.new %>% group_by(gene) %>% summarize(n=n()) %>% filter(n==4) %>% as.data.frame()

# ## //
# a3 <- encc.case2.melt %>% group_by(gene, severity) %>% summarize(zero_count=sum(iszero), total_count=n_distinct(cell)) %>%
# 	as.data.frame()
# a3$missing.rate <- a3$zero_count / a3$total_count

# a3.new <- a3 %>% filter(missing.rate <0.5)
# a3.new.in <- a3.new %>% group_by(gene) %>% summarize(n=n()) %>% filter(n==4) %>% as.data.frame()


# ## // compare mRNA expression of included genes between severity groups

# encc.case1.melt.in <- encc.case1.melt %>% filter(gene %in% a2.new.in$gene)
# encc.case1.gene.mean <- encc.case1.melt.in %>% group_by(gene, severity) %>% summarize(exp_mean = mean(exp))
# encc.case1.gene.wide <- reshape2::dcast(encc.case1.gene.mean, gene~severity, value.var="exp_mean")

# encc.case2.melt.in <- encc.case2.melt %>% filter(gene %in% a3.new.in$gene)
# encc.case2.gene.mean <- encc.case2.melt.in %>% group_by(gene, severity) %>% summarize(exp_mean = mean(exp))
# encc.case2.gene.wide <- reshape2::dcast(encc.case2.gene.mean, gene~severity, value.var="exp_mean")

# ##
# tmp <- encc.case1.gene.wide
# rownames(tmp) <- tmp$gene
# tmp <- tmp[,2:5]

# pdf(file.path(pub.data.dir, "encc.case1.gene.mrna.heatmap.pdf"))
# Heatmap(tmp, cluster_columns = FALSE)
# dev.off()

# ##
# tmp <- encc.case2.gene.wide
# rownames(tmp) <- tmp$gene
# tmp <- tmp[,2:5]

# pdf(file.path(pub.data.dir, "encc.case2.gene.mrna.heatmap.pdf"))
# Heatmap(tmp, cluster_columns = FALSE)
# dev.off()

# ##



# ## ...............
# ## // pca of single encc cells

# encc.data.case1.in <- encc.data[encc.case1.gene.wide$gene,]
# encc.data.case2.in <- encc.data[encc.case2.gene.wide$gene,]

# tmp <- t(as.matrix((encc.data.case1.in)))
# encc.data.case1.pca <- PCA(tmp,  graph = FALSE)
# scree.plot <- fviz_screeplot(encc.data.case1.pca, addlabels = TRUE, ncp = ncol(tmp))

# pdf(file.path(pub.data.dir, "encc.case1.pca.scree.plot.pdf"))
# print(scree.plot)
# dev.off()

# ## 

# tmp <- t(as.matrix((encc.data.case2.in)))
# res.pca <- PCA(tmp,  graph = FALSE)
# scree.plot <- fviz_screeplot(res.pca, addlabels = TRUE, ncp = ncol(tmp))

# pdf(file.path(pub.data.dir, "encc.case2.pca.scree.plot.pdf"))
# print(scree.plot)
# dev.off()


## -------------------
## 




