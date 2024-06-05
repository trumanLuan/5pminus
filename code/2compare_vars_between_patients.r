## genetics of Cri du chat and HSCR

require(readr)
require(dplyr)
require(stringr)
require(Seurat)
# library(SummarizedExperiment) 

require(foreach)
require(doParallel) ## for parallel computation

proj.dir <- "/home/luanyz0418/project/clinical_cases/maojiao"
data.dir <- file.path(proj.dir, 'data')
pub.data.dir <- file.path(proj.dir, 'public_data')

hscr.gene <- readxl::read_xlsx(path=file.path(proj.dir,"HSCR-genes.xlsx"), sheet="Gene-based") %>% as.data.frame()


## --------------------------
## load data
## full table of snv, cnv and sv of case1 and case2
## --------------------------

## data description

## case1 snvs, rare frequency, exonic regions, filter out snvs present in parents
## case1 cnvs, filter out overlap regions with dgvMerged and cnvs present in parents, re-annotation to refGenes and cytoband.

## case2 snvs, rare frequency, exonic regions
## case2 cnvs, filter out overlap regions with dgvMerged, re-annotation to refGenes and cytoband
## case2 svs, filter out overlap regions with dgvMerged, re-annotation to refGenes, and cytoband

snv.case1 <- read.table(file.path(data.dir, 'case1/case1.snv.final.tsv'), header=T, sep='\t', as.is=T )
cnv.case1 <- read.table(file.path(data.dir, 'case1/case1.cnv.final.tsv'), header=T, sep='\t', as.is=T )

snv.case2 <- readr::read_rds( file.path(data.dir, "case2/case2.snvs.filter.rds") ) %>% as.data.frame()
cnv.case2 <- readr::read_rds( file.path(data.dir, "case2/case2.cnvs.filter.rds") )
sv.case2 <- readr::read_rds( file.path(data.dir, "case2/case2.svs.filter.rds") )

##

snv.case2.new <- snv.case2[,c(1:20, 148:149)]
snv.case2.new <- snv.case2.new %>% filter(`Func.refGene` %in% c('exonic', 'exonic;splicing', 'splicing', 'UTR3', 'UTR5', 'UTR5;UTR3')) %>% 
						filter(`ExonicFunc.refGene` %in% c("frameshift insertion", "nonsynonymous SNV", "stopgain", "stoploss")) %>% filter(Coverage >20)

sv.case2.new <- sv.case2 %>% filter(gene_region %in% c("exonic", "splicing", "UTR3", "UTR5", "UTR5;UTR3")) # 25548 rows
sv.case2.new$SU <- sapply(sv.case2.new$info, function(x) stringr::str_split(x, ":")[[1]][2])
sv.case2.new$SU <- as.numeric(sv.case2.new$SU)
sv.case2.new <- sv.case2.new %>% filter(SU > 20) # 264 events

cnv.case2.new <- cnv.case2 %>% filter(gene_region %in% c("exonic", "splicing", "UTR3", "UTR5", "UTR5;UTR3")) # 692 rows

# d.cnv.sub <- d.cnvs %>% filter(gene_region %in% c("exonic", "splicing", "UTR3", "UTR5")) # 692 rows
# d.sv.sub <- d.svs %>% filter(gene_region %in% c("exonic", "splicing", "UTR3", "UTR5", "UTR5;UTR3")) # 25548 rows


## --------------------------
## direct comparison of snvs, cnvs and svs between case1 and case2.
## using full tables
## --------------------------

## ..............
##  directly compare snvs

tmp <- stringr::str_replace_all(snv.case2$CHROM, 'chr', '')
snv.case2.list <- stringr::str_c(tmp, ':', snv.case2$POS, ':', snv.case2$REF,':',snv.case2$ALT)
snv.case1.list <- snv.case1$Otherinfo1

tmp <- Reduce(intersect, list(snv.case1.list, snv.case2.list))
snv.case1[snv.case1$Otherinfo1 %in% tmp,]
snv.case2[snv.case2.list %in% tmp, 1:10]

## ..............
##  directly compare cnv regions

cnv.case2.new <- rbind(cnv.case2[,c('chr', 'start', 'end', 'svtype', 'svlen', 'gene_region', 'genes', 'cytoBand')],
	sv.case2[,c('chr', 'start', 'end', 'svtype', 'svlen', 'gene_region', 'genes', 'cytoBand')])

cnv.case1.loss <- cnv.case1 %>% filter(svtype.final == 'DEL')
cnv.case1.gain <- cnv.case1 %>% filter(svtype.final == 'DUP')

cnv.case2.loss <- cnv.case2.new %>% filter(svtype %in% c('loss', 'DEL'))
cnv.case2.gain <- cnv.case2.new %>% filter(svtype %in% c('gain', "DUP"))

cm.cnv.loss <- bedtoolsr::bt.intersect(a= cnv.case1.loss[,1:3], b= cnv.case2.loss[,1:3])
cm.cnv.gain <- bedtoolsr::bt.intersect(a= cnv.case1.gain[,1:3], b= cnv.case2.gain[,1:3])

cm.cnv.loss <- unique(cm.cnv.loss)

## annotate to refGenes and cytoband
cm.cnv.loss$V4 = 0
cm.cnv.loss$V5 = 0
write.table(cm.cnv.loss, file.path(data.dir, "common.cnv.loss.between.case1.and.case2.tsv"), quote=F, sep='\t', row.names=F, col.names=F)

## // shell code
## // annovar annotation of sv events 
cd /home/luanyz0418/utility/annovar
anv_input=/home/luanyz0418/project/clinical_cases/maojiao/data/common.cnv.loss.between.case1.and.case2.tsv
./annotate_variation.pl -build hg19 -out common.cnv $anv_input humandb/    ## annotate to refGene.
./annotate_variation.pl -regionanno -build hg19 -out common.cnv -dbtype cytoBand $anv_input humandb/   ## annotate to cytoBand


## reorganize cm.cnv data frame by adding refgene and cytoband info
tmp.cytoband <- read.table("/home/luanyz0418/utility/annovar/common.cnv.hg19_cytoBand", header=F, sep='\t', as.is=T)
tmp.refgene <- read.table("/home/luanyz0418/utility/annovar/common.cnv.variant_function", header=F, sep='\t', as.is=T)

colnames(cm.cnv.loss) <- c("chr", 'start', 'end', 'svlen', 'sv.region')
cm.cnv.loss$svlen <- cm.cnv.loss$end - cm.cnv.loss$start
cm.cnv.loss$`sv.region` <- tmp.refgene[,1]
cm.cnv.loss$genes <- tmp.refgene[,2]
cm.cnv.loss$cytoband <- tmp.cytoband[,2]

cm.cnv.loss$chr.arm <- ifelse(grepl('p', cm.cnv.loss$cytoband), 'p', 'q')

cm.cnv.loss %>% group_by(chr, chr.arm) %>% summarize(sum.svlen = sum(svlen))

# check protein-coding genes located in 5p 
a <- cm.cnv.loss %>% filter(chr == '5' & chr.arm == 'p')
a.protein <- a %>% filter( sv.region %in% c('exonic', 'intronic', 'UTR3', 'UTR5'))
a.ncrna <-  a %>% filter( sv.region %in% c('ncRNA_exonic', 'ncRNA_intronic' ))

## check proein-coding genes in all the shared cnv loss regions
b.protein <- cm.cnv.loss %>% filter(sv.region %in% c('exonic', 'intronic', 'UTR3', 'UTR5'))
b.ncrna <- cm.cnv.loss %>% filter( sv.region %in% c('ncRNA_exonic', 'ncRNA_intronic' ))

b.protein.genes <- unique( unlist(sapply(b.protein$genes, function(x) stringr::str_split(x, ',')[[1]]) ) )
b.protein.genes <- b.protein.genes[! grepl("lins0\\)", b.protein.genes)]
tmp <- unique( sapply(b.protein.genes, function(x) stringr::str_split(x, '\\(')[[1]][1]) )

table(tmp %in% hscr.gene$Gene)
