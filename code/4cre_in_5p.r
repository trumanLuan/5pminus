
library(reshape2)
library(dplyr)

	library(ggVennDiagram)
	library(ggplot2)

library(ComplexHeatmap)
library(circlize)

library(karyoploteR)
library(ChIPpeakAnno)
library(GenomicFeatures)

## ------------------------------------------
proj.dir <- "/home/luanyz0418/project/clinical_cases/maojiao/public_data/scATAC"
data.dir <- file.path(proj.dir, "data/")
genome.dir <- "/home/luanyz0418/project/clinical_cases/maojiao/genome_data"

hscr.gene <- readxl::read_xlsx(path=file.path("/home/luanyz0418/project/clinical_cases/maojiao","HSCR-genes.xlsx"), sheet="Gene-based") %>% as.data.frame()

## ------------------------------------------
## total number of CRE in human genome
cre.total <- read.table(file.path(data.dir, 'adult/CRE_total/cCRE_hg38.tsv'), as.is=T, header=F, sep='\t')
colnames(cre.total) <- c('chr', 'start', 'end', 'class', 'present_in_fetal', "present_in_adult", "module")

## cytoband info for hg38
cytoband.info <- readr::read_tsv(file.path(data.dir, 'adult/CRE_total/cytoBand.txt'), col_names=F)
colnames(cytoband.info) <- c('chr', 'start', 'end', 'cytoband', 'desc')
cytoband.info.subset <- cytoband.info %>% filter(chr == "chr5") %>% as.data.frame() ## end of 5p 48,800,000

cytoband.hg19 <- readr::read_tsv("/home/luanyz0418/project/cffDNA/genome_files/hg19/cytoBand_ucsc_hg19.txt", col_names=F)
colnames(cytoband.hg19) <- c('chr', 'start', 'end', 'cytoband', 'desc')
cytoband.hg19.subset <- cytoband.hg19 %>% filter(chr == "chr5") %>% as.data.frame() ## end of 5p 48,400,000

## check cre list located in 5p
cre.subset <- subset(cre.total, chr == 'chr5' )
cre.subset <- cre.subset[cre.subset$end < 48800000,]
cre.subset$index <- stringr::str_c(cre.subset$chr, ":", cre.subset$start, ":", cre.subset$end)

## stat of class
table(cre.subset$class)

	## stat of present_in_fetal tissues

	x <- list(fetal=subset(cre.subset, present_in_fetal =='yes')$index, 
	          adult=subset(cre.subset, present_in_adult == "yes")$index )

	ggVennDiagram(x, label_alpha = 0, set_color = 'black') + 
	  scale_fill_distiller(palette = "white")
	ggsave(file.path(proj.dir, "Venn_plot_of_cre_overlap_between_adult_and_fetal.pdf"))



## ------------------------------------------
## cell type specitif usage of CREs, figure 

cre.list <-  read.table(file.path(data.dir, 'adult/CRE_by_cell_type/CRE.bed'), as.is=T, sep='\t', header=F)
colnames(cre.list) <- c('chr', 'start', 'end')
cre.list$index <- stringr::str_c(cre.list$chr, ":", cre.list$start, ':', cre.list$end)

## input the matrix of CRE and cell types
cre.mtx <- read.table(file.path(data.dir, 'adult/CRE_by_cell_type/matrix.tsv'), skip=3, as.is=T, sep=' ', header=F)
colnames(cre.mtx) <- c("cre", 'celltype')
cre.mtx$value <- rep(1, nrow(cre.mtx))

focus.cre.id <- which( cre.list$index %in% cre.subset$index )
cre.mtx.subset <- cre.mtx %>% filter(cre %in% focus.cre.id)

cre.mtx.new <-  dcast(cre.mtx.subset, cre ~ celltype, value.var="value")
cre.mtx.new <- cre.mtx.new[, 2:ncol(cre.mtx.new)]

rownames(cre.mtx.new) <- cre.subset$index

cre.mtx.new[which(is.na(cre.mtx.new), arr.ind = T)] <- 0

a <- rowSums(cre.mtx.new)
summary(a)


	## heatmap of CRE usage
	pdf(file.path(proj.dir, "heatmap_5p_cres_in_all_celltypes.pdf"), height=8, width=8)
	Heatmap(cre.mtx.new, show_row_names = F, show_column_names = F, col=colorRamp2(c(0, 1), c("white", "red")) )
	dev.off()


## --------------------------------------------
## focus on the intestinal cells

cell.index <- c(4, 5, 18, 53, 66, 67, 94, 96, 103:111, 144, 151, 154, 212)
names(cell.index) <- c("T_CD8", 'T_CD4', 'Fibro_GI', "Enteric_Neuron", 'Plasma_B', 'Memory_B', 'Sm_Ms_Colon', 'Sm_Ms_GI',
                       'Colon_Epithelial_1', 'Enterocyte', 'Colon_Goblet', 'SI_Goblet', 'Colon_Epithelial_2', 'Colon_Epithelial_3',
                       "Enterochromaffin", "Tuft", 'Paneth', "Fetal_enteroendocrine", 'Fetal_Enteric_neuron', 'Fetal_ENteric_Glia',
                       "Fetal_Fibro_GI")
cre.mtx.GI <- cre.mtx.new[,cell.index]
cre.mtx.GI <- cre.mtx.GI[rowSums(cre.mtx.GI)>0,] #3849 CREs
colnames(cre.mtx.GI) <- names(cell.index)

	## heatmap of CRE usage in intestinal cells.
		pdf(file.path(proj.dir, "heatmap_5p_cres_in_gut_celltypes.pdf"), height=8, width=8)
		Heatmap(cre.mtx.GI, show_row_names = F, show_column_names = T, col=colorRamp2(c(0, 1), c("white", "red")) )
		dev.off()

## CREs in "Enteric_Neuron", 'Fetal_Enteric_neuron', 'Fetal_ENteric_Glia',

test.mtx <- cre.mtx.GI[,c("Enteric_Neuron", 'Fetal_Enteric_neuron', 'Fetal_ENteric_Glia')]
table( rowSums(test.mtx) > 0 )

## -------------------------------------------
## analyze the CRE-gene links in enteric neurons

## adult tissues

tmp.dat <- read.table(file.path(data.dir, 'adult/enteric_neuron.tsv'), as.is=T, sep='\t', header=T)
tmp.dat <- tmp.dat[grepl('chr5:', tmp.dat$cCRE),]


tmp.dat$chr <- sapply(tmp.dat$cCRE, function(x) stringr::str_split(x, ":")[[1]][1])
tmp.dat$chr_start <- sapply(tmp.dat$cCRE, function(x) {
  a = stringr::str_split(x, ":")[[1]][2]
  stringr::str_split(a, '-')[[1]][1]
  })
tmp.dat.subset <- subset(tmp.dat, chr == 'chr5'  )
tmp.dat.subset <- subset(tmp.dat.subset, as.numeric(chr_start) < 48800000 )

tmp.dat.subset$gene <- sapply(tmp.dat.subset$`Gene.Name`, function(x) stringr::str_split(x, ':')[[1]][1])
tmp.dat.subset$gene_ens <- sapply(tmp.dat.subset$`Gene.Name`, function(x) stringr::str_split(x, ':')[[1]][2])

gene.ens.short <- sapply(tmp.dat.subset$gene_ens, function(x) stringr::str_split(x, '\\.')[[1]][1])


		## load hg38 gtf file, and check gene start and end.
		library(GenomicFeatures)
		
		gtf.tx <- GenomicFeatures::makeTxDbFromGFF(file=file.path(genome.dir, 'hg38.ensGene.gtf.gz'), format = 'gtf')
		exon.by.tr <- GenomicFeatures::exonsBy(gtf.tx, by='gene', use.name=T) 
		cds.by.tr <- GenomicFeatures::cdsBy(gtf.tx, by='tx',use.name=T)
		cds.by.tr <- cds.by.tr[names(cds.by.tr) %in% unique(gene.ens.short)]

		cds.combined <- NULL
		for(i in 1:length(cds.by.tr)){
			tmp <- cds.by.tr[[i]] %>% as.data.frame()
			tmp$geneid <- rep(names(cds.by.tr)[[i]], nrow(tmp))
			cds.combined <- rbind(cds.combined, tmp)
		}

		all(cds.combined$start < 48800000)
		all(cds.combined$end < 48800000)

		## functional enrichment analysis of target genes
		library(gprofiler2)

		enrich.res <- gost(unique(tmp.dat.subset$gene), organism = "hsapiens")


# tmp.dat.subset$index <- tmp.dat.subset$cCRE
# tmp.dat.subset$index <- stringr::str_replace_all(tmp.dat.subset$index, "-", ":")
# table( tmp.dat.subset$index %in% rownames(cre.mtx.GI) )


FINAL.SUBSET <- subset(tmp.dat.subset, grepl("GDNF", Gene.Name))
write.table(FINAL.SUBSET, "E:/FOCUS_in_5p.txt",quote=F,sep='\t',row.names = F, col.names = T)


		## plot links 
		# library(karyoploteR)
		# start.regs <- toGRanges(data.frame("chr1", 20e6, 30e6))
		# end.regs <- toGRanges(data.frame("chr3", 50e6, 55e6))

		# kp <- plotKaryotype()
		# kpPlotLinks(kp, data=start.regs, data2=end.regs)

			start.df.chr <- sapply(tmp.dat.subset$cCRE, function(x) stringr::str_split(x, ":")[[1]][1] )
			start.df.reg <- sapply(tmp.dat.subset$cCRE, function(x) stringr::str_split(x, ":")[[1]][2] )
			start.df.reg2 <- sapply(start.df.reg, function(x) stringr::str_split(x, '-')[[1]])
			start.df <- data.frame(chr=start.df.chr, start = as.numeric( t(start.df.reg2)[,1] ), end = as.numeric( t(start.df.reg2)[,2] ))

			end.df.chr <- sapply(tmp.dat.subset$Promoter, function(x) stringr::str_split(x, ":")[[1]][1] )
			end.df.reg <- sapply(tmp.dat.subset$Promoter, function(x) stringr::str_split(x, ":")[[1]][2] )
			end.df.reg2 <- sapply(end.df.reg, function(x) stringr::str_split(x, '-')[[1]])
			end.df <- data.frame(chr=end.df.chr, start = as.numeric( t(end.df.reg2)[,1] ), end = as.numeric( t(end.df.reg2)[,2] ))

			start.regs <- toGRanges(start.df)
			end.regs <- toGRanges(end.df)

			pdf(file.path(proj.dir, 'link_map_adult_cre_and_genes3.pdf'), width=140, height=8)
			kp <- plotKaryotype(genome='hg19', chromosomes = 'chr5', plot.type = 2)
			kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")
			kpAddCytobandLabels(kp)

			kpDataBackground(kp, r0=0, r1=0.15)
			# kpPlotRegions(kp, data=exon.by.tr[['ENSG00000168621']])
			cnv.case1 <- read.table(file.path("/home/luanyz0418/project/clinical_cases/maojiao/data", 'case1/case1.cnv.final.tsv'), header=T, sep='\t', as.is=T )
			cnv.regions.case1 <- cnv.case1 %>% filter(chr.final == 5) 
			cnv.regions.case1 <- cnv.regions.case1[,c('chr.final', 'start.final', 'end.final')]
			cnv.regions.case1$chr.final <- stringr::str_c("chr",cnv.regions.case1$chr.final)
			kpPlotRegions(kp, data=toGRanges(data.frame(chr=cnv.regions.case1$chr.final, start=cnv.regions.case1$start.final, end=cnv.regions.case1$end.final)), r0=0, r1=0.15) # for patient 1

			kpDataBackground(kp, r0=0.2, r1=0.35)
			cnv.case2 <- readr::read_rds( file.path("/home/luanyz0418/project/clinical_cases/maojiao/data", "case2/case2.cnvs.filter.rds") )
			cnv.case2.regions <- cnv.case2 %>% filter(svtype == 'loss' & chr == 5)
			cnv.case2.regions$chr <- stringr::str_c("chr",cnv.case2.regions$chr)

			sv.case2 <- readr::read_rds( file.path("/home/luanyz0418/project/clinical_cases/maojiao/data", "case2/case2.svs.filter.rds") )
			sv.case2$SU <- sapply(sv.case2$info, function(x) stringr::str_split(x, ":")[[1]][2])
			sv.case2$SU <- as.numeric(sv.case2$SU)
			sv.case2.regions <- sv.case2 %>% filter(SU >20 & svtype == "DEL" & chr == 5) # 57 reginos

			kpPlotRegions(kp, data=toGRanges(data.frame(chr=cnv.case2.regions$chr, start=cnv.case2.regions$start, end=cnv.case2.regions$end)), r0=0.2, r1=0.35) # for patient 2
			kpPlotLinks(kp, data=start.regs, data2=end.regs, arch.height = 0.8, r0=0.4, r1=1)
			dev.off()
## ///////////

## fetal cells

# coacc.score <- readr::read_csv(file.path(data.dir, "fetal/GSE149683_File_S5.Cicero_coaccessibility_scores_by_cell_type.csv"), col_names=T)
# sub.rows <- grepl("chr5-", coacc.score$Peak1) | grepl("chr5-", coacc.score$Peak2)
# coacc.score <- coacc.score[sub.rows,] 
# coacc.score <- coacc.score %>% as.data.frame()

# coacc.score <- coacc.score[,c("Peak1", 'Peak2', "ENS glia_intestine", "ENS neurons_intestine")]

# sub.rows2 <- apply(coacc.score[,3:4], 1, function(x) all(is.na(x)))
# coacc.score <- coacc.score[!sub.rows2,]
# readr::write_rds(coacc.score, file.path(data.dir, "fetal/coaccessibility_score_chr5_subset.rds"))

coacc.score <- readr::read_rds( file.path(data.dir, "fetal/coaccessibility_score_chr5_subset.rds") )

colnames(coacc.score) <- c('peak1', 'peak2', 'glia', 'neuron')

coacc.score$peak1_chr <- sapply(coacc.score$peak1, function(x) stringr::str_split(x, '-')[[1]][1])
coacc.score$peak1_start <- sapply(coacc.score$peak1, function(x) stringr::str_split(x, '-')[[1]][2])
coacc.score$peak1_end <- sapply(coacc.score$peak1, function(x) stringr::str_split(x, '-')[[1]][3])

coacc.score$peak2_chr <- sapply(coacc.score$peak2, function(x) stringr::str_split(x, '-')[[1]][1])
coacc.score$peak2_start <- sapply(coacc.score$peak2, function(x) stringr::str_split(x, '-')[[1]][2])
coacc.score$peak2_end <- sapply(coacc.score$peak2, function(x) stringr::str_split(x, '-')[[1]][3])


coacc.score <- coacc.score %>% filter(as.numeric(peak1_start) < 48400000) %>% filter(as.numeric(peak2_start) < 48400000)
dup.rows <- t( apply(coacc.score[,1:2], 1, sort) )
dup.rows <- as.data.frame(dup.rows)
dup.rows.id <- duplicated(dup.rows)
coacc.score <- coacc.score[!dup.rows.id,]
		
	coacc.score$d12 <- as.numeric(coacc.score$peak1_end) - as.numeric(coacc.score$peak2_start)
	coacc.score$d21 <- as.numeric(coacc.score$peak2_end) - as.numeric(coacc.score$peak1_start)
	coacc.score$d <- apply(coacc.score[,c('d12','d21')], 1, function(x){
		if(abs(x[1]) > abs(x[2])) tmp.d = abs(x[2])
		if(abs(x[2]) > abs(x[1])) tmp.d = abs(x[1])
		tmp.d
		})

	summary(coacc.score$d) ## distance between CREs and target genes


## anntate the intervals to genes and promoters
library(ChIPpeakAnno)

df <- coacc.score[,c("peak1_chr", "peak1_start", "peak1_end")]
colnames(df) <- c('chr', "start", 'end')
peak.list.peak1 <- makeGRangesFromDataFrame(df)

df <- coacc.score[,c("peak2_chr", "peak2_start", "peak2_end")]
colnames(df) <- c('chr', "start", 'end')
peak.list.peak2 <- makeGRangesFromDataFrame(df)

peak.annot <- GenomicFeatures::makeTxDbFromGFF(file=file.path(genome.dir, 'hg19.refGene.gtf.gz'), format = 'gtf')
seqlevels(peak.annot) <- 'chr5' ## set active chrs in current TxDb object
y <- GenomicFeatures::transcriptsBy(peak.annot, by='gene')


peak1.annot <- annotatePeakInBatch(myPeakList = peak.list.peak1, AnnotationData = unlist(y)) 
peak2.annot <- annotatePeakInBatch(myPeakList = peak.list.peak2, AnnotationData = unlist(y)) 

names(peak1.annot) <- NULL; peak1.annot <- peak1.annot %>% as.data.frame()
names(peak2.annot) <- NULL; peak2.annot<- peak2.annot %>% as.data.frame()
peak.annot <- rbind(peak1.annot, peak2.annot)
peak.annot$index <- stringr::str_c(peak.annot$seqnames, '-', peak.annot$start, '-', peak.annot$end)

peak1.position <- sapply(coacc.score$peak1, function(x){
	tmp <- subset(peak.annot, index %in% x)
	stringr::str_c(unique(tmp$insideFeature), collapse=';')
	})
peak2.position <- sapply(coacc.score$peak2, function(x){
	tmp <- subset(peak.annot, index %in% x)
	stringr::str_c(unique(tmp$insideFeature), collapse=';')
	})
peak1.gene <- sapply(coacc.score$peak1, function(x){
	tmp <- subset(peak.annot, index %in% x)
	stringr::str_c(unique(tmp$feature), collapse=';')
	})
peak2.gene <- sapply(coacc.score$peak2, function(x){
	tmp <- subset(peak.annot, index %in% x)
	stringr::str_c(unique(tmp$feature), collapse=';')
	})

coacc.score$peak1_position <- peak1.position
coacc.score$peak2_position <- peak2.position
coacc.score$peak1_gene <- peak1.gene
coacc.score$peak2_gene <- peak2.gene

coacc.score$pos_pair <- stringr::str_c(coacc.score$peak1_position, '_', coacc.score$peak2_position)

coacc.score.inuse <- coacc.score %>% filter(! pos_pair %in% c("downstream_downstream", "downstream_upstream", "upstream_downstream", "upstream_upstream"))
coacc.score.inuse <- coacc.score.inuse %>% filter(as.numeric(peak1_start) < 48400000) %>% filter(as.numeric(peak2_start) < 48400000)


	## functional enrichement 
	library(gprofiler2)

	enrich.res <- gost(unique(c(coacc.score.inuse$peak1_gene, coacc.score.inuse$peak2_gene)), organism = "hsapiens")
	dim(enrich.res$result)


## check whether HSCR-geens were involved
table(coacc.score.inuse$peak1_gene %in% hscr.gene$Gene)
table(coacc.score.inuse$peak2_gene %in% hscr.gene$Gene)

tmp.sub <- coacc.score.inuse$peak1_gene %in% hscr.gene$Gene | coacc.score.inuse$peak2_gene %in% hscr.gene$Gene
coacc.score.inuse[tmp.sub,]



		## plot links of fetal data
			fetal.start.df.chr <- sapply(coacc.score$peak1, function(x) stringr::str_split(x, "-")[[1]][1] )
			fetal.start.df.reg <- sapply(coacc.score$peak1, function(x) stringr::str_split(x, "-")[[1]][2] )
			fetal.start.df.reg2 <- sapply(coacc.score$peak1, function(x) stringr::str_split(x, '-')[[1]][3] )
			fetal.start.df <- data.frame(chr=fetal.start.df.chr, start = as.numeric( (fetal.start.df.reg) ), end = as.numeric( (fetal.start.df.reg2) ))

			fetal.end.df.chr <- sapply(coacc.score$peak2, function(x) stringr::str_split(x, "-")[[1]][1] )
			fetal.end.df.reg <- sapply(coacc.score$peak2, function(x) stringr::str_split(x, "-")[[1]][2] )
			fetal.end.df.reg2 <- sapply(coacc.score$peak2, function(x) stringr::str_split(x, '-')[[1]][3])
			fetal.end.df <- data.frame(chr=fetal.end.df.chr, start = as.numeric( (fetal.end.df.reg) ), end = as.numeric( (fetal.end.df.reg2) ))

			fetal.start.regs <- toGRanges(fetal.start.df)
			fetal.end.regs <- toGRanges(fetal.end.df)

			pdf(file.path(proj.dir, 'link_map_fetal_cre_and_genes.pdf'), width=140, height=8)
			kp <- plotKaryotype(genome='hg19', chromosomes = 'chr5', plot.type = 2)
			kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1,minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")
			kpAddCytobandLabels(kp)
			kpPlotLinks(kp, data=fetal.start.regs, data2=fetal.end.regs, arch.height=0.8)
			dev.off()