
## conda env: base

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

library(foreach)
library(monocle)

library(clusterProfiler) # functional enrichment anlaysis
library(msigdbr)
library(enrichplot)

# library(factoextra)
# library("FactoMineR")

library(ComplexHeatmap)
library(circlize)

hscr.gene <- readxl::read_xlsx(path=file.path(proj.dir,"HSCR-genes.xlsx"), sheet="Gene-based") %>% as.data.frame()


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
table(gene.case1 %in% rownames( gutcell.data) ) # 72 out of 99  in
table(gene.case2 %in% rownames( gutcell.data ) ) # 508 out of 732  in
table(gene.shared %in% rownames( gutcell.data ) ) # 7 out of 8


## differential analysis between cell types
Idents(gutcell.obj) <- 'Integrated_05'

cell.use <- c('Adult Glia', "Branch A1 (iMN)", "Branch A2 (IPAN/IN)", "Branch A3 (IPAN/IN)", "Branch A4 (IN)",
			"Branch B1 (eMN)", "Branch B2 (eMN)", "Branch B3 (IPAN)", "Differentiating glia", "ENCC/glia Progenitor",
			"Glia 1 (DHH+)", "Glia 2 (ELN+)", "Glia 3 (BCAN+)", "cycling ENCC/glia")


# res.combine <- NULL
# for(x in cell.use ) {
# 	cat(date(), x, '\n')
# 	res.marker <- FindMarkers(gutcell.obj, ident.1 = x, ident.2 = NULL, only.pos = FALSE)
# 	res.marker$cluster <- rep(x, nrow(res.marker))
# 	res.marker$gene <- rownames(res.marker)
# 	res.combine <- rbind(res.combine, res.marker)
# }
# 
# 
# readr::write_rds(res.combine, file.path(pub.data.dir, "gutcellalas_neuron_findmarkers.rds"))

res.combine <- readr::read_rds( file.path(pub.data.dir, "gutcellalas_neuron_findmarkers.rds") )

### ```` only check in positive markers

res.combine.markers <- res.combine %>% filter(avg_log2FC >1 & p_val_adj < 0.05 & pct.1 >0.3)

res.combine.markers.wide <- dcast(res.combine.markers, gene ~ cluster, value.var = 'avg_log2FC')

test.wide <- res.combine.markers.wide[,2:ncol(res.combine.markers.wide)]
rownames(test.wide) <- res.combine.markers.wide$gene
test.wide[which(is.na(test.wide), arr.ind =T)] <- 0

pdf(file.path(data.dir, '../test.pdf') ) ## fig 2A

row.anno <- rowAnnotation(patient1=ifelse(rownames(test.wide) %in% gene.case1, 1, 0), 
	patient2=ifelse(rownames(test.wide) %in% gene.case2, 1, 0) )

ComplexHeatmap::Heatmap(test.wide, col = colorRamp2(c(0, 6), c("white", "deeppink")), 
	border_gp = gpar(col = "black", lty = 2), show_row_names=FALSE, right_annotation = row.anno  )
dev.off()

a1 = sort( gene.case1[gene.case1 %in% rownames(test.wide)] )

a2 = sort(gene.case2[gene.case2 %in% rownames(test.wide)] )

##

ego <- enrichGO(gene=res.combine.markers.wide$gene, universe=rownames(gutcell.data), OrgDb='org.Hs.eg.db', keyType="SYMBOL", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, readable=T)
ego.df <- as.data.frame(ego)
write.table(ego.df, file.path(data.dir, '../GO_terms_enteric_neuronal_cell_genes.tsv'),quote=F, sep='\t', row.names=F, col.names=T)

##
mdf = msigdbr(species="Homo sapiens")

c8_t2g <- mdf %>% filter(gs_cat == "C8") %>% 
  dplyr::select(gs_name, gene_symbol)

enrich.c8 <- enricher(gene=res.combine.markers.wide$gene, universe=rownames(gutcell.data), TERM2GENE=c8_t2g)
enrich.c8.df = as.data.frame(enrich.c8)
write.table(enrich.c8.df, file.path(data.dir, '../enricher_celltypesignatures_enteric_neuronal_cell_genes.tsv'),quote=F, sep='\t', row.names=F, col.names=T)

enrich.r <- enrichplot::pairwise_termsim(enrich.c8)

p1 <- enrichplot::treeplot(enrich.r, showCategory=20)

ggsave(file.path(data.dir, "../treeplot_enrichr_celltypesignatures_encc.pdf"), width=15, height=8) ## figure 2B

##
res.combine.markers.case1 <- res.combine.markers %>% filter(gene %in% gene.case1)
res.combine.markers.case2 <- res.combine.markers %>% filter(gene %in% gene.case2)

### ```end

# res.combine.sig <- res.combine %>% filter(p_val_adj < 0.05)
# res.combine.sig$fc_direction <- ifelse(res.combine.sig$avg_log2FC > 0, 'up', 'down')

# 	a <- res.combine.sig %>% group_by(cluster, fc_direction) %>% summarize(n=n()) %>% as.data.frame()
# 	a$fc_direction <- factor(a$fc_direction, levels=c('up', 'down'))

# 	p <- ggplot(a, aes(x=cluster, y=n, fill=fc_direction)) + 
# 	  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) +
# 	  theme_light() + coord_flip()
# 	ggsave(file.path(pub.data.dir, "plot_encc_degs_in_gutcellatlas.pdf"))


# 	## 
# 	table(gene.case1 %in% res.combine.sig$gene) # 58 in 99 genes, differentially expressed
# 	table(gene.case2 %in% res.combine.sig$gene) # 381 in 732 genes, differentially expressed

# 	res.combine.sig.sub1 <- res.combine.sig %>% filter(gene %in% c(gene.case1))
# 	a <- res.combine.sig.sub1 %>% group_by(cluster, fc_direction) %>% summarize(n=n()) %>% as.data.frame()
# 	a$fc_direction <- factor(a$fc_direction, levels=c('up', 'down'))

# 	p <- ggplot(a, aes(x=cluster, y=n, fill=fc_direction)) + 
# 	  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) +
# 	  theme_light() + coord_flip()
# 	ggsave(file.path(pub.data.dir, "plot_encc_degs_in_gutcellatlas_case1.pdf"))

# 	##
# 	res.combine.sig.sub2 <- res.combine.sig %>% filter(gene %in% c(gene.case2))
# 	a<-res.combine.sig.sub2 %>% group_by(cluster, fc_direction) %>% summarize(n=n()) %>% as.data.frame()
# 	a$fc_direction <- factor(a$fc_direction, levels=c('up', 'down'))

# 	p <- ggplot(a, aes(x=cluster, y=n, fill=fc_direction)) + 
# 	  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) +
# 	  theme_light() + coord_flip()
# 	ggsave(file.path(pub.data.dir, "plot_encc_degs_in_gutcellatlas_case2.pdf"))

## 


## ----------------------------
## single cell pseudotime ordering of gut neurons
## monocle 2
## --------------------------

## subset single cells from Seurat object
Idents(gutcell.obj) <- "Integrated_05"
gutcell.obj.sub <- subset(gutcell.obj, idents = cell.use)

	## t-SNE plot 
	pdf(file.path(data.dir, "pca.seurat.gutcell.obj.neurons.pdf"))
	g <- DimPlot(object = gutcell.obj.sub, label=TRUE, reduction = 'umap')
	print(g)
	dev.off()


	pdf(file.path(data.dir, "pca.seurat.gutcell.obj.neurons.tsne.pdf"))
	g <- DimPlot(object = gutcell.obj.sub, label=TRUE, reduction = 'tsne')
	print(g)
	dev.off()


## //
## using all neuron cells to create new CellDataSet object

pd <- new("AnnotatedDataFrame", data = gutcell.obj.sub@meta.data)
fd <- new("AnnotatedDataFrame", data = gutcell.obj.sub@assays$RNA@meta.features)
cds <- newCellDataSet(gutcell.obj.sub@assays$RNA@counts, phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())


# # (Required) Estimate size factors and dispersions， ～ 10 mins
# cds <- estimateSizeFactors(cds)
# cds <- estimateDispersions(cds)

# ## cluster single cells
# # cds <- clusterCells(cds)

# ## order single cells in pseudotime along a trajectory
# disp_table <- dispersionTable(cds)
# ordering_genes <- subset(disp_table, mean_expression >= 0.1)
# cds <- setOrderingFilter(cds, ordering_genes)
# cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree') # ~ time consuming, around 4 hrs.
# cds <- orderCells(cds) # ~ 20 mins, 

# readr::write_rds(cds, file.path(data.dir, "monocle_obj_after_orderCells.rds"))
cds <- readr::read_rds( file.path(data.dir, "monocle_obj_after_orderCells.rds") )


	## plot single cell trajectory
	pdf(file.path(data.dir, "plot.cell.trajectory.by.celltype.pdf"),width=10, height=10)
	g <- plot_cell_trajectory(cds, color_by = "Integrated_05")
	print(g)
	dev.off()

	pdf(file.path(data.dir, "plot.cell.trajectory.by.celltype.split.pdf"),width=15, height=15)
	g <- plot_cell_trajectory(cds, color_by = "Integrated_05") + facet_wrap(~`Integrated_05`, nrow=3)
	print(g)
	dev.off()

	## 
	pdf(file.path(data.dir, "plot.cell.trajectory.by.agegroup.pdf"),width=10, height=10)
	g <- plot_cell_trajectory(cds, color_by = "Age_group")
	print(g)
	dev.off()

	pdf(file.path(data.dir, "plot.cell.trajectory.by.agegroup.split.pdf"),width=10, height=10)
	g <- plot_cell_trajectory(cds, color_by = "Age_group") + facet_wrap(~ `Age_group`, nrow=2)
	print(g)
	dev.off()

	## 
	pdf(file.path(data.dir, "plot.cell.trajectory.by.Pseudotime.pdf"),width=10, height=10)
	g <- plot_cell_trajectory(cds, color_by = "Pseudotime")
	print(g)
	dev.off()

		## 
	pdf(file.path(data.dir, "plot.cell.trajectory.by.state.pdf"),width=10, height=10)
	g <- plot_cell_trajectory(cds, color_by = "State")
	print(g)
	dev.off()

## // select genes shared between two cases and single cell dataset
select.genes.case1 <- Reduce(intersect, list(gene.case1, fData(cds)$gene_ids )) # 72 shared genes.
select.genes.case2 <- Reduce(intersect, list(gene.case2, fData(cds)$gene_ids )) # 508 shared genes.


# ## Finding genes that changes as a function of pseudotime
# diff.genes.pseudotime.case1 <- differentialGeneTest(cds[select.genes.case1,], fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 1, verbose = 1) # ~1 min
# sig.genes.pseudotime.case1 <- subset(diff.genes.pseudotime.case1, qval < 0.05) #  50 genes.

# diff.genes.pseudotime.case2 <- differentialGeneTest(cds[select.genes.case2,], fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 1, verbose = 1) # , ~
# sig.genes.pseudotime.case2 <- subset(diff.genes.pseudotime.case2, qval < 0.05) # 302 genes.


# ## finding genes that distinguish cell type or state
# diff.genes.celltype.case1 <- differentialGeneTest(cds[select.genes.case1,], fullModelFormulaStr = "~`Integrated_05`", cores = 1, verbose = 1) # ~
# sig.genes.celltype.case1 <- subset(diff.genes.celltype.case1, qval < 0.05) # 58

# diff.genes.celltype.case2 <- differentialGeneTest(cds[select.genes.case2,], fullModelFormulaStr = "~`Integrated_05`", cores = 1, verbose = 1) # ~
# sig.genes.celltype.case2 <- subset(diff.genes.celltype.case2, qval < 0.05) # 377


## analyze genes with changes along the branches of cell trajectory
# beam.methods <- c('duplicate', 'expression', 'cluster')
# beam.res.case1 <- lapply(beam.methods, function(method) {
# 	BEAM(cds[select.genes.case1,], branch_point = 1, progenitor_method = method, cores=3)
# })

# monocle::BEAM(cds[select.genes.case1,], branch_point=1, cores=1, progenitor_method = 'duplicate')

fData(cds)$gene_short_name <- fData(cds)$gene_ids


## plot_genes_branched_heatmap
a.case1 = monocle::BEAM(cds[select.genes.case1,], branch_point=1, cores=1, progenitor_method = 'duplicate') 
a.case1 <- a.case1[order(a.case1$qval),]
a.case1 <- a.case1[,c('gene_short_name', 'pval', 'qval')]

	pdf(file.path(pub.data.dir, "plot.genes.branched.heatmap.case1.pdf"),width=10, height=10)
	plot_genes_branched_heatmap(cds[row.names(subset(a.case1, qval < 0.05)),],
                                          branch_point = 1,
                                          num_clusters = 4,
                                          # add_annotation_row= cbind(gene.snv.annot.case1, gene.sv.annot.case1),
                                          # add_annotation_row= row.annot,
                                          cores = 1,
                                          use_gene_short_name = T,
                                          show_rownames = T,
                                          return_heatmap = F)
	dev.off()


# ttt = cds[row.names(subset(a.case1, qval < 0.05)),]
## a.case1 genes in res.combine.markers and hscr.genes
a.case1 %>% filter(qval < 0.05) %>% filter( (gene_short_name %in% hscr.genes) & (gene_short_name %in% res.combine.markers$gene) )

# row.annot <- data.frame(var_info = gene.sv.annot.case1$var_info)


a.case2 = monocle::BEAM(cds[select.genes.case2,], branch_point=1, cores=1, progenitor_method = 'duplicate') 
a.case2 <- a.case2[order(a.case2$qval),]
a.case2 <- a.case2[,c('gene_short_name', 'pval', 'qval')]

	pdf(file.path(pub.data.dir, "plot.genes.branched.heatmap.case2.pdf"),width=10, height=15)
	plot_genes_branched_heatmap(cds[row.names(subset(a.case2, qval < 0.05)),],
                                          branch_point = 1, hclust_method = "average",
                                          num_clusters = 5,
                                          cores = 1,
                                          use_gene_short_name = T, 
                                          show_rownames = T,
                                          return_heatmap = F)
	dev.off()

##
hscr.genes <- Reduce(intersect, list(rownames(cds), hscr.gene$Gene))
	pdf(file.path(pub.data.dir, "plot.genes.branched.heatmap.hscr.genes.pdf"),width=10, height=15)
	plot_genes_branched_heatmap(cds[hscr.genes,],
                                          branch_point = 1, hclust_method = "complete",
                                          num_clusters = 6,
                                          cores = 1,
                                          use_gene_short_name = T, 
                                          show_rownames = T,
                                          return_heatmap = F)
	dev.off()

	## annotate genes in the plot_genes_branches_heatmap with mutation type information
	gene.snv.annot.case1 <- foreach(x = row.names(subset(a.case1, qval < 0.05)), .combine='rbind' ) %do% {
			tmp.snv.info <- snv.case1.select[grepl(x, snv.case1.select$`Gene.refGene`),]
			# tmp.cnv.info <- cnv.case1.select[grepl(x, cnv.case1.select$`genes.final`),]

			if(nrow(tmp.snv.info) > 0){
				var_type = 'SNV'
				var_region = stringr::str_c(unique(tmp.snv.info$`Func.refGene`), collapse = ',')
				var_info = stringr::str_c(unique(tmp.snv.info$`ExonicFunc.refGene`), collapse = ',')
			}else{
				var_type = 'n.a.'
				var_region = 'n.a.'
				var_info = 'n.a.'
			}

			data.frame(var_type = var_type, var_region = var_region, var_info = var_info)

		}

	gene.sv.annot.case1 <- foreach(x = row.names(subset(a.case1, qval < 0.05)), .combine='rbind' ) %do% {
			tmp.sv.info <- cnv.case1.select[grepl(x, cnv.case1.select$`genes.final`),]
			# tmp.cnv.info <- cnv.case1.select[grepl(x, cnv.case1.select$`genes.final`),]

			if(nrow(tmp.sv.info) > 0){
				var_type = 'SV'
				var_region = stringr::str_c(unique(tmp.sv.info$`region.final`), collapse = ',')
				var_info = stringr::str_c(unique(tmp.sv.info$`svtype.final`), collapse = ',')
			}else{
				var_type = 'n.a.'
				var_region = 'n.a.'
				var_info = 'n.a.'
			}

			data.frame(var_type = var_type, var_region = var_region, var_info = var_info)

		}


## highlight some genes, and plot_genes_branched_pseudotime
	pdf(file.path(pub.data.dir, "plot.genes.branched.pseudotime.known.hscr.genes.pdf"),width=10, height=15)
	plot_genes_branched_pseudotime(cds[c("RET", 'NRG1', "ERBB2", "ERBB3", "ITGB1"),],
                       branch_point = 1,
                       color_by = "Pseudotime",
                       ncol = 2)
	dev.off()


	pdf(file.path(pub.data.dir, "plot.genes.branched.pseudotime.by.state.known.hscr.genes.pdf"),width=10, height=15)
	plot_genes_branched_pseudotime(cds[c("RET", 'NRG1', "ERBB2", "ERBB3", "ITGB1"),],
                       branch_point = 1,
                       color_by = "State",
                       ncol = 2)
	dev.off()

	pdf(file.path(pub.data.dir, "plot.genes.branched.pseudotime.by.celltype.known.hscr.genes.pdf"),width=10, height=15)
	plot_genes_branched_pseudotime(cds[c("RET", 'NRG1', "ERBB2", "ERBB3", "ITGB1"),],
                       branch_point = 1,
                       color_by = "Integrated_05",
                       ncol = 2)
	dev.off()


##
    pdf(file.path(pub.data.dir, "plot.genes.branched.pseudotime.by.state.case1.pdf"),width=10, height=10)
	plot_genes_branched_pseudotime(cds[c("CDH10", "CDH12", "NPR3", "DAG1"),],
                       branch_point = 1,
                       color_by = "State",
                       ncol = 2)
	dev.off()



## --------------
## load ENCC single cell rnaseq data
## --------------

hscr.data = file.path(pub.data.dir, "HSCR.seurat.Rdata")
load(hscr.data) #HSCR.seurat
rm(hscr.data)

encc.obj <- HSCR.seurat

encc.meta <- encc.obj@meta.data
encc.meta$cell <- rownames(encc.meta)

encc.data <- encc.obj@assays$RNA@data

## check overlap between gene vars and genes included in single cell data
table(gene.case1 %in% rownames( encc.data) )
table(gene.case2 %in% rownames( encc.data ) )
table(gene.shared %in% rownames( encc.data ) ) # 7 out of 8


## differntial analysis between severity groups
Idents(encc.obj) <- 'severity'
res.marker <- FindAllMarkers(encc.obj, only.pos = FALSE)

table(gene.case1 %in% res.marker$genes)
table(gene.case2 %in% res.marker$genes)








