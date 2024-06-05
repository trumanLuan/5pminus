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

## ---------------
## load hscr genes.

hscr.gene <- readxl::read_xlsx(path=file.path(proj.dir,"HSCR-genes.xlsx"), sheet="Gene-based") %>% as.data.frame()

## -----------------------
## load case1 sample SNP, indel, cnv and sv 
## -----------------------


## load case1 data, old data from Zuo-xiaoyu
case1.d.variants <- readr::read_tsv(file = file.path(data.dir, 'case1/WGS.82.HMX.deeply_filtered.hg19_multianno.txt'))
case1.d.var.vcf <- read.table(file = file.path(data.dir, 'case1/WGS.82.HMX.deeply_filtered.renameID.vcf'), header=F)
colnames(case1.d.var.vcf) <- c('CHROM', 'POS', "ID", "REF", "ALT", 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'HMX', 'HMX_F', 'HMX_M')
case1.d.var <- merge(case1.d.variants, case1.d.var.vcf, by.x = "Otherinfo1", by.y='ID')

case1.d.var.sub <- case1.d.var %>% select(Otherinfo1, Chr, Start, End, Ref, Alt, cytoBand, Func.refGene, Gene.refGene, GeneDetail.refGene, ExonicFunc.refGene, AAChange.refGene, HMX, HMX_F, HMX_M)
case1.d.var.sub$HMX <- stringr::str_replace_all(case1.d.var.sub$HMX, "\\|", "/")
case1.d.var.sub$HMX_F <- stringr::str_replace_all(case1.d.var.sub$HMX_F, "\\|", "/")
case1.d.var.sub$HMX_M <- stringr::str_replace_all(case1.d.var.sub$HMX_M, "\\|", "/")

case1.d.var.sub$HMX <- sapply(case1.d.var.sub$HMX, function(x) stringr::str_split(x, ":")[[1]][1])
case1.d.var.sub$HMX_F <- sapply(case1.d.var.sub$HMX_F, function(x) stringr::str_split(x, ":")[[1]][1])
case1.d.var.sub$HMX_M <- sapply(case1.d.var.sub$HMX_M, function(x) stringr::str_split(x, ":")[[1]][1])

case1.d.var.sub <- case1.d.var.sub[case1.d.var$`1000g2015aug_all` < 0.01,]
case1.d.var.sub <- case1.d.var.sub[!is.na(case1.d.var.sub$HMX),]

is.out <- case1.d.var.sub$HMX == case1.d.var.sub$HMX_F | case1.d.var.sub$HMX == case1.d.var.sub$HMX_M
case1.d.snv.final <- case1.d.var.sub[!is.out,]

is.out <- case1.d.snv.final$HMX == '0/1' & case1.d.snv.final$HMX_F == '1/1'
case1.d.snv.final <- case1.d.snv.final[! is.out,]

is.out <- case1.d.snv.final$HMX == './.'
case1.d.snv.final <- case1.d.snv.final[! is.out,]

is.out <- case1.d.snv.final$HMX == '0/1' & case1.d.snv.final$HMX_M == '1/1'
case1.d.snv.final <- case1.d.snv.final[! is.out,]

write.table(case1.d.snv.final,file.path(data.dir, 'case1/case1.snv.final.tsv'), quote=F, sep='\t', row.names=F, col.names=T)


## ------------------------
## load case1 cnv data
case1.d.cnv <- readxl::read_xlsx(path = file.path(data.dir, 'case1/clinSV/results/SV-CNV.RARE_PASS_GENE.xlsx'))
case1.d.cnv$idindex <- stringr::str_c(case1.d.cnv$ID,':',case1.d.cnv$pedInfo1)

cnv.pos <- case1.d.cnv %>% select(idindex, LOCATION) %>% as.data.frame()
tmp.chr <- sapply(cnv.pos$LOCATION, function(x) stringr::str_split(x, ':')[[1]][1])
tmp.pos <- sapply(cnv.pos$LOCATION, function(x) stringr::str_split(x, ':')[[1]][2])
tmp.start <- sapply(tmp.pos, function(x) stringr::str_split(x, '-')[[1]][1])
tmp.end <- sapply(tmp.pos, function(x) stringr::str_split(x, '-')[[1]][2])

cnv.annovar.input <- data.frame(chr=tmp.chr, start=tmp.start, end=tmp.end, ref=rep(0, length(tmp.chr)), alt=rep(0, length(tmp.chr)))
cnv.annovar.input <- unique(cnv.annovar.input)
# write.table(cnv.annovar.input, file.path(data.dir, 'case1/annovar.input.cnv.results.tsv'), quote=F,sep='\t',row.names=F,col.names=F)

# ## // shell bash
# ## // annovar annotation of cnv events
# cd /home/luanyz0418/utility/annovar
# anv_input=/home/luanyz0418/project/clinical_cases/maojiao/data/case1/annovar.input.cnv.results.tsv
# ./annotate_variation.pl -regionanno -build hg19 -out maojiao_case1 -dbtype cytoBand $anv_input humandb/

cytoband.anno <- read.table(file.path(data.dir, 'case1/maojiao_case1.hg19_cytoBand'),header=F, sep='\t', as.is=T)
colnames(cytoband.anno) <- c('dbtype', 'cytoband', 'chr', 'start', 'end', 'ref', 'alt')
cytoband.anno$id <- stringr::str_c(cytoband.anno$chr, ':', cytoband.anno$start, '-', cytoband.anno$end)

case1.d.cnv <- as.data.frame(case1.d.cnv)
case1.d.cnv$cytoband <- sapply(case1.d.cnv$LOCATION, function(x) {
	x.out = cytoband.anno[cytoband.anno$id %in% x, 'cytoband']
	if(length(x.out)==0 ) x.out='n.a.'
	x.out
}
	)

write.table(case1.d.cnv, file.path(data.dir, 'case1/cnv.results.annovar.combined.tsv'), quote=F,sep='\t',row.names=F,col.names=T)


# ## 
# tmp <- case1.d.cnv %>% filter(SAMPLE == 'HMX') %>% select(LOCATION, SVTYPE)
# tmp$chr <- sapply(tmp$LOCATION, function(x) stringr::str_split(x, ":")[[1]][1])
# tmp.loca <- sapply(tmp$LOCATION, function(x) stringr::str_split(x, ":")[[1]][2])
# tmp$start <- sapply(tmp.loca, function(x) stringr::str_split(x, "-")[[1]][1])
# tmp$end <- sapply(tmp.loca, function(x) stringr::str_split(x, "-")[[1]][2])
# tmp <- tmp %>% select(chr, start, end, SVTYPE)
# tmp <- tmp[!grepl(',',tmp$start),]
# write.table(tmp, file.path(data.dir, 'case1/cnv.results.HMX.bed'), quote=F,sep='\t',row.names=F,col.names=F)

# tmp <- case1.d.cnv %>% filter(SAMPLE == 'HMX-F') %>% select(LOCATION, SVTYPE)
# tmp$chr <- sapply(tmp$LOCATION, function(x) stringr::str_split(x, ":")[[1]][1])
# tmp.loca <- sapply(tmp$LOCATION, function(x) stringr::str_split(x, ":")[[1]][2])
# tmp$start <- sapply(tmp.loca, function(x) stringr::str_split(x, "-")[[1]][1])
# tmp$end <- sapply(tmp.loca, function(x) stringr::str_split(x, "-")[[1]][2])
# tmp <- tmp %>% select(chr, start, end, SVTYPE)
# tmp <- tmp[!grepl(',',tmp$start),]
# write.table(tmp, file.path(data.dir, 'case1/cnv.results.HMX_F.bed'), quote=F,sep='\t',row.names=F,col.names=F)

# tmp <- case1.d.cnv %>% filter(SAMPLE == 'HMX-M') %>% select(LOCATION, SVTYPE)
# tmp$chr <- sapply(tmp$LOCATION, function(x) stringr::str_split(x, ":")[[1]][1])
# tmp.loca <- sapply(tmp$LOCATION, function(x) stringr::str_split(x, ":")[[1]][2])
# tmp$start <- sapply(tmp.loca, function(x) stringr::str_split(x, "-")[[1]][1])
# tmp$end <- sapply(tmp.loca, function(x) stringr::str_split(x, "-")[[1]][2])
# tmp <- tmp %>% select(chr, start, end, SVTYPE)
# tmp <- tmp[!grepl(',',tmp$start),]
# write.table(tmp, file.path(data.dir, 'case1/cnv.results.HMX_M.bed'), quote=F,sep='\t',row.names=F,col.names=F)


## // check cnv overlapping regions between daughter and parents
# Reduce(intersect, list(1:12,12:20) )
cnv.hmx <- readr::read_tsv(file.path(data.dir, 'case1/cnv.results.HMX.bed'), col_names=F)
cnv.hmxf <- readr::read_tsv(file.path(data.dir, 'case1/cnv.results.HMX_F.bed'), col_names=F)
cnv.hmxm <- readr::read_tsv(file.path(data.dir, 'case1/cnv.results.HMX_M.bed'), col_names=F)

col.name <- c('chr', 'start', 'end', 'svtype')
colnames(cnv.hmx) <- col.name
colnames(cnv.hmxf) <- col.name
colnames(cnv.hmxm) <- col.name


mark_shared_regions <- function(data1, data2){
	cl <- makeCluster(20)
	registerDoParallel(cl)

	data1.renew <- foreach(ri = 1:nrow(data1), .combine=rbind) %dopar% {
	#for(ri in 1:nrow(data1)){
		cat(ri, '\n')

		curr.chr <- subset(data2, chr %in% data1$chr[ri])

		test.overlap <- c()

		for(j in 1:nrow(curr.chr)){
			if(data1$svtype[ri] == curr.chr$svtype[j]){
				curr.overlap <- Reduce(intersect, list(data1$start[ri]:data1$end[ri], curr.chr$start[j]:curr.chr$end[j]))
				if(length(curr.overlap) > 0){
					test.overlap <- c(test.overlap, TRUE)
				}else{
					test.overlap <- c(test.overlap, FALSE)
				}
			}else{
				test.overlap <- c(test.overlap,FALSE)
			}
		}

		cbind(data1[ri,],overlap=any(test.overlap))
	}

	stopCluster(cl)
	data1.renew
}

cnv.hmx.overlap.hmxf <- mark_shared_regions(cnv.hmx, cnv.hmxf)
cnv.hmx.overlap.hmxm <- mark_shared_regions(cnv.hmx, cnv.hmxm)
is.overlap <- cnv.hmx.overlap.hmxf$overlap | cnv.hmx.overlap.hmxm$overlap

tmp <- cnv.hmx.overlap.hmxf[is.overlap,] %>% filter(overlap) %>% select(chr, start, end, svtype)
cnv.hmx.rm.father = bedtoolsr::bt.subtract(a= tmp, b = cnv.hmxf)
cnv.hmx.rm.parent = bedtoolsr::bt.subtract(a= cnv.hmx.rm.father, b= as.data.frame(cnv.hmxm) )
colnames(cnv.hmx.rm.parent) <- c('chr', 'start', 'end', 'svtype')

cnv.hmx.uniq.final <- rbind(cnv.hmx.rm.parent, cnv.hmx[!is.overlap,])


## // remove cnv overlapping regions between case1 and dgv database

dgv.cnvs <- readr::read_tsv('/home/luanyz0418/utility/annovar/humandb/hg19_dgvMerged.txt', col_names=F)
dgv.cnvs.subcols <- dgv.cnvs[,c(2:4,11)]
colnames(dgv.cnvs.subcols) <- c('chr', 'start', 'end', 'svtype')
dgv.cnvs.subcols$chr <- stringr::str_replace_all(dgv.cnvs.subcols$chr, 'chr', '')

##
dgv.cnvs.del <- dgv.cnvs.subcols %>% filter(svtype %in% c('deletion', 'loss', 'gain+loss'))
dgv.cnvs.dup <- dgv.cnvs.subcols %>% filter(svtype %in% c("duplication", "gain", "gain+loss"))


##
cnv.hmx.uniq.deletion <- cnv.hmx.uniq.final %>% filter(svtype == 'DEL')
cnv.hmx.uniq.duplication <- cnv.hmx.uniq.final %>% filter(svtype == "DUP")

cnv.hmx.uniq.del.rm.dgv <- bedtoolsr::bt.subtract(a= cnv.hmx.uniq.deletion, b= dgv.cnvs.del)
cnv.hmx.uniq.dup.rm.dgv <- bedtoolsr::bt.subtract(a= cnv.hmx.uniq.duplication, b= dgv.cnvs.dup)
colnames(cnv.hmx.uniq.del.rm.dgv) <- c('chr', 'start', 'end', 'svtype')
colnames(cnv.hmx.uniq.dup.rm.dgv) <- c('chr', 'start', 'end', 'svtype')


## ......combine
cnv.hmx.uniq.final.renew <- rbind(cnv.hmx.uniq.del.rm.dgv, cnv.hmx.uniq.dup.rm.dgv, subset(cnv.hmx.uniq.final, svtype == "BND"))
cnv.hmx.uniq.final.renew$len <- cnv.hmx.uniq.final.renew$end - cnv.hmx.uniq.final.renew$start

write.table(cnv.hmx.uniq.final.renew, file.path(data.dir, "case1/cnv.hmx.uniq.final.tsv"),quote=F,sep='\t', row.names=F, col.names=T)

cnv.hmx.uniq.final.renew$svtype <- 0
cnv.hmx.uniq.final.renew$len <- 0
write.table(cnv.hmx.uniq.final.renew, file.path(data.dir, "case1/cnv.hmx.uniq.final.annovar.input.tsv"),quote=F,sep='\t', row.names=F, col.names=F)



# ## // shell bash
# ## // annovar annotation of cnv events, according to 
# cd /home/luanyz0418/utility/annovar
# anv_input=/home/luanyz0418/project/clinical_cases/maojiao/data/case1/cnv.hmx.uniq.final.annovar.input.tsv
# ./annotate_variation.pl -build hg19 -out maojiao_case1 $anv_input humandb/

## check erbb3 regions
tt <- read.table( file.path(data.dir, "case1/cnv.hmx.uniq.final.tsv"), header=T, sep='\t', as.is=T)
tt2 <- read.table( file.path(data.dir, "case1/maojiao_case1.cnv.hmx.uniq.final.variant_function"), header=F, sep='\t', as.is=T)
colnames(tt2) <- c('region', 'genes', 'chr', 'start', 'end', 'v6', 'v7')
case1.d.cnv.final.uniq <- cbind(tt2[,c('chr', 'start', 'end','region', 'genes')], tt[,c('svtype', 'len')])

case1.d.cnv.subcols <- case1.d.cnv %>% filter(pedInfo1 == "daughter") %>% select(LOCATION, GT, DRF, DRA, SVTYPE) %>% as.data.frame()
case1.d.cnv.subcols$chr <- sapply(case1.d.cnv.subcols$LOCATION, function(x) stringr::str_split(x, ":")[[1]][1])
tmp <- sapply(case1.d.cnv.subcols$LOCATION, function(x) stringr::str_split(x, ":")[[1]][2])
case1.d.cnv.subcols$start <- sapply(tmp, function(x) stringr::str_split(x, "-")[[1]][1])
case1.d.cnv.subcols$end <- sapply(tmp, function(x) stringr::str_split(x, "-")[[1]][2])

case1.d.cnv.subcols <- case1.d.cnv.subcols %>% select(chr, start, end, GT, DRF, DRA, SVTYPE) %>% filter(SVTYPE %in% c("DEL", "DUP"))

test <- bedtoolsr::bt.intersect(a= case1.d.cnv.final.uniq, b= case1.d.cnv.subcols, wa=T, wb=T)
colnames(test) <- c('chr.final', 'start.final', 'end.final', 'region.final', 'genes.final', 'svtype.final', 'svlen.final', 
	'chr.hmx.raw', 'start.hmx.raw', 'end.hmx.raw', 'GT.hmx.raw', 'DRF.hmx.raw', 'DRA.hmx.raw', 'svtype.hmx.raw')

write.table(test, file.path(data.dir, "case1/hmx.cnv.final.table2.tsv"), quote=F,sep='\t',row.names=F, col.names=T)

## map the renew regions of case1 cnvs to cytoband

case1.cnv.final <- read.table(file.path(data.dir, "case1/hmx.cnv.final.table2.tsv"), header=T, sep='\t', as.is=T)
test <- case1.cnv.final[,1:5]
test$`region.final` <- 0
test$`genes.final` <- 0
write.table(test, file.path(data.dir, "case1/hmx.cnv.final.table2.simple.annovar.input.bed"), quote=F,sep='\t',row.names=F, col.names=F)


# ## // shell bash
# ## // annovar annotation of cnv events
# cd /home/luanyz0418/utility/annovar
# anv_input=/home/luanyz0418/project/clinical_cases/maojiao/data/case1/hmx.cnv.final.table2.simple.annovar.input.bed
# ./annotate_variation.pl -regionanno -build hg19 -out maojiao_case1.hmx.cnv.final -dbtype cytoBand $anv_input humandb/
# mv maojiao_case1.* ~//project/clinical_cases/maojiao/data/case1

cytoband.anno <- read.table(file.path(data.dir, 'case1/maojiao_case1.hmx.cnv.final.hg19_cytoBand'),header=F, sep='\t', as.is=T)
case1.cnv.final$cytoband <- cytoband.anno[,2]
write.table(case1.cnv.final, file.path(data.dir, 'case1/case1.cnv.final.tsv'), quote=F, sep='\t',row.names=F, col.names=T)


## // check the genes located in the exonic CNVs.
case1.cnv.genes.exonic <- case1.cnv.final %>% filter(region.final %in% c('exonic', 'splicing'))
case1.snv.genes.exonic <- read.table(file.path(data.dir, 'case1/case1.snv.final.tsv'), header=T, sep='\t', as.is=T)

cnv.gene.list <- unique(case1.cnv.genes.exonic$`genes.final`)
cnv.gene.list2 <- unique(unlist(sapply(cnv.gene.list, function(x) unlist(stringr::str_split(x, ',')[[1]]) ) ) )

snv.gene.list <- unique(case1.snv.genes.exonic$`Gene.refGene`)
snv.gene.list2 <- unique(unlist(sapply(snv.gene.list, function(x) unlist(stringr::str_split(x, ';')[[1]]) ) ) )


table(snv.gene.list2 %in% hscr.gene$Gene )
table(cnv.gene.list2 %in% hscr.gene$Gene )


## -----------------------
## load case2 cnv and sv data
## -----------------------


d.snv <- readr::read_tsv(file = file.path(data.dir, 'case2/analysis_report/4-Annotation/XZD/XZD.annotation.xls'))
d.cnv <- readr::read_tsv(file = file.path(data.dir, 'case2/analysis_report/5-CNV/XZD/XZD.CNVs.annotation.xls'))
d.sv <- readr::read_tsv(file = file.path(data.dir, 'case2/analysis_report/6-SV/XZD/XZD.SVs.annotation.xls'))

# ## ......................

# ## annotate d.snv with `1000g2015aug_all`
# d.snv.tmp<- d.snv[,c('CHROM', 'POS', 'POS', 'REF', 'ALT')]
# write.table(d.snv.tmp, file.path(data.dir, 'case2/case2.snvs.all.annovar.input.tsv'), quote=F,sep='\t', row.names=F, col.names=F)


# ## -- shell bash
# cd /home/luanyz0418/utility/annovar
# anv_input=/home/luanyz0418/project/clinical_cases/maojiao/data/case2/case2.snvs.all.annovar.input.tsv
# ./annotate_variation.pl -filter -build hg19 -out maojiao_case2.snv.all -dbtype 1000g2015aug_all -maf 0.01 $anv_input humandb/ 

# ## -- filter and select cols.

# annovar.filter <- read.table('/home/luanyz0418/utility/annovar/maojiao_case2.snv.all.hg19_ALL.sites.2015_08_filtered', header=F, sep='\t', as.is=T)

# d.snv$index <- stringr::str_c(d.snv$CHROM, ':', d.snv$POS, ':', d.snv$REF, ':', d.snv$ALT)
# annovar.filter$index <- stringr::str_c(annovar.filter[,1], ':', annovar.filter[,2], ':', annovar.filter[,4], ':', annovar.filter[,5])


# # d.variants.new <- d.variants[,c(1:19, 138:148)] %>% as.data.frame()

# ## check the annotated genes of all genomic variants
# d.snv.sub <- d.snv %>% filter(index %in% annovar.filter$index) %>%
# 	filter(Func.refGene %in% c('exonic', 
# 		'exonic;splicing', 
# 		'ncRNA_exonic',
# 		'ncRNA_exonic;splicing',
# 		'ncRNA_splicing',
# 		'ncRNA_UTR5',
# 		'splicing',
# 		'UTR3', 
# 		'UTR5', 
# 		'UTR5;UTR3')) 

# write_rds(d.snv.sub, file.path(data.dir, "case2/case2.snvs.filter.rds"))


## ......................
## annotate and compare cnvs to dgvMerged database

## load dgvMerged cnvs
dgv.cnvs <- readr::read_tsv('/home/luanyz0418/utility/annovar/humandb/hg19_dgvMerged.txt', col_names=F)
dgv.cnvs.subcols <- dgv.cnvs[,c(2:4,11)]
colnames(dgv.cnvs.subcols) <- c('chr', 'start', 'end', 'svtype')
dgv.cnvs.subcols$chr <- stringr::str_replace_all(dgv.cnvs.subcols$chr, 'chr', '')

##
dgv.cnvs.del <- dgv.cnvs.subcols %>% filter(svtype %in% c('deletion', 'loss', 'gain+loss'))
dgv.cnvs.dup <- dgv.cnvs.subcols %>% filter(svtype %in% c("duplication", "gain", "gain+loss"))


## load case2 cnvs.

## included cnvs mapped to exonic regions
d.cnv.sub <- d.cnv %>% filter(Func.refGene %in% c('exonic', 'ncRNA_exonic', 'UTR5')) %>% as.data.frame()

d.cnv.del <- d.cnv.sub %>% filter(Alteration == 'loss') %>% select(Chromosome, Start, End, PredictedCopyNumber, Alteration)
d.cnv.dup <- d.cnv.sub %>% filter(Alteration == "gain") %>% select(Chromosome, Start, End, PredictedCopyNumber, Alteration)

## subtract cnv regions shared with dgvMerged

d.cnv.del.rm.dgv <- bedtoolsr::bt.subtract(a= d.cnv.del, b= dgv.cnvs.del)
d.cnv.dup.rm.dgv <- bedtoolsr::bt.subtract(a= d.cnv.dup, b= dgv.cnvs.dup)

colnames(d.cnv.del.rm.dgv) <- c('chr', 'start', 'end', 'cnv_nubmer', 'svtype')
colnames(d.cnv.dup.rm.dgv) <- c('chr', 'start', 'end', 'cnv_nubmer', 'svtype')

d.cnv.del.rm.dgv$svlen <- d.cnv.del.rm.dgv$end - d.cnv.del.rm.dgv$start
d.cnv.dup.rm.dgv$svlen <- d.cnv.dup.rm.dgv$end - d.cnv.dup.rm.dgv$start

## annotate to refGene with annovar
tmp<- rbind(d.cnv.del.rm.dgv[,1:5], d.cnv.dup.rm.dgv[,1:5] )
tmp[,4]= 0
tmp[,5] = 0
write.table(tmp, file.path(data.dir, 'case2/case2.cnvs.all.annovar.input.tsv'), quote=F, sep='\t', row.names=F, col.names=F)

## // shell code
## // annovar annotation of cnv events using --geneanno mode 
cd /home/luanyz0418/utility/annovar
anv_input=/home/luanyz0418/project/clinical_cases/maojiao/data/case2/case2.cnvs.all.annovar.input.tsv
./annotate_variation.pl -build hg19 -out maojiao_case2.cnvs $anv_input humandb/
./annotate_variation.pl -regionanno -build hg19 -out maojiao_case2.cnvs -dbtype cytoBand $anv_input humandb/ ## annotate to cytoBand


## combine columns
d.cnv.annot.refgene <- read.table('/home/luanyz0418/utility/annovar/maojiao_case2.cnvs.variant_function', header=F, sep='\t', as.is=T)
d.cnv.annot.cytoband <- read.table('/home/luanyz0418/utility/annovar/maojiao_case2.cnvs.hg19_cytoBand', header=F, sep='\t', as.is=T)

d.cnv.rm.dgv <- rbind(d.cnv.del.rm.dgv, d.cnv.dup.rm.dgv)
d.cnv.rm.dgv$gene_region <- d.cnv.annot.refgene[,1]
d.cnv.rm.dgv$genes <- d.cnv.annot.refgene[,2]
d.cnv.rm.dgv$cytoBand <- d.cnv.annot.cytoband[,2]


write_rds(d.cnv.rm.dgv, file.path(data.dir, "case2/case2.cnvs.filter.rds"))

## ......................
## annotate and compare SVs to dgvMerged database

d.sv <- d.sv %>% as.data.frame()

d.sv$svtype <- sapply(d.sv$INFO, function(x) stringr::str_split(x, ";")[[1]][1])
d.sv$svtype <- stringr::str_replace_all(d.sv$svtype, "SVTYPE=", "")
d.sv$svlen <- sapply(d.sv$INFO, function(x) stringr::str_split(x, ";")[[1]][3])
d.sv$svlen <- stringr::str_replace_all(d.sv$svlen, "SVLEN=", "")

d.sv.sub <- d.sv %>% filter(svtype %in% c("DEL", "DUP"))
d.sv.sub$svlen <- as.numeric(d.sv.sub$svlen)

d.sv.sub.new <- data.frame(chr = d.sv.sub[,1], start=d.sv.sub$POS, end=d.sv.sub$POS+abs(d.sv.sub$svlen)-1, 
	ref=d.sv.sub$REF, svtype=d.sv.sub$svtype, svlen=d.sv.sub$svlen, format=d.sv.sub$FORMAT, info=d.sv.sub$XZD)


## remove shared regions with dgvMerged database

d.sv.del <- d.sv.sub.new %>% filter(svtype == 'DEL') 
d.sv.dup <- d.sv.sub.new %>% filter(svtype == "DUP") 

d.sv.del$chr <- stringr::str_replace_all(d.sv.del$chr, 'chr', '')
d.sv.dup$chr <- stringr::str_replace_all(d.sv.dup$chr, 'chr', '')

d.sv.del.rm.dgv <- bedtoolsr::bt.subtract(a= d.sv.del, b= dgv.cnvs.del)
d.sv.dup.rm.dgv <- bedtoolsr::bt.subtract(a= d.sv.dup, b= dgv.cnvs.dup)

colnames(d.sv.del.rm.dgv) <- c('chr', 'start', 'end', 'ref', 'svtype', 'svlen', 'format', 'info')
colnames(d.sv.dup.rm.dgv) <- c('chr', 'start', 'end', 'ref', 'svtype', 'svlen', 'format', 'info')

d.sv.del.rm.dgv$svlen <- d.sv.del.rm.dgv$end - d.sv.del.rm.dgv$start + 1
d.sv.dup.rm.dgv$svlen <- d.sv.dup.rm.dgv$end - d.sv.dup.rm.dgv$start + 1


## annotate sv regions to refgene and cytoband
tmp <- rbind(d.sv.del.rm.dgv, d.sv.dup.rm.dgv)[,1:5]
tmp[,4] = 0
tmp[,5] = 0

write.table(tmp, file.path(data.dir, 'case2/case2.svs.all.annovar.input.tsv'), quote=F, sep='\t', row.names=F, col.names=F)

## // shell code
## // annovar annotation of sv events 
cd /home/luanyz0418/utility/annovar
anv_input=/home/luanyz0418/project/clinical_cases/maojiao/data/case2/case2.svs.all.annovar.input.tsv
./annotate_variation.pl -build hg19 -out maojiao_case2.svs $anv_input humandb/    ## annotate to refGene.
./annotate_variation.pl -regionanno -build hg19 -out maojiao_case2.svs -dbtype cytoBand $anv_input humandb/   ## annotate to cytoBand

## combine tables

d.sv.annot.refgene <- read.table('/home/luanyz0418/utility/annovar/maojiao_case2.svs.variant_function', header=F, sep='\t', as.is=T)
d.sv.annot.cytoband <- read.table('/home/luanyz0418/utility/annovar/maojiao_case2.svs.hg19_cytoBand', header=F, sep='\t', as.is=T)

d.sv.rm.dgv <- rbind(d.sv.del.rm.dgv, d.sv.dup.rm.dgv)
d.sv.rm.dgv$gene_region <- d.sv.annot.refgene[,1]
d.sv.rm.dgv$genes <- d.sv.annot.refgene[,2]
d.sv.rm.dgv$cytoBand <- d.sv.annot.cytoband[,2]


write_rds(d.sv.rm.dgv, file.path(data.dir, "case2/case2.svs.filter.rds"))

## ''''''''''''''''''''
## ''''''''''''''''''''


d.snvs <- readr::read_rds( file.path(data.dir, "case2/case2.snvs.filter.rds") ) %>% as.data.frame()
d.cnvs <- readr::read_rds( file.path(data.dir, "case2/case2.cnvs.filter.rds") )
d.svs <- readr::read_rds( file.path(data.dir, "case2/case2.svs.filter.rds") )


## statistics of exonic hits
table(d.snvs$`Func.refGene`) 
table(d.cnvs$gene_region)
table(d.svs$gene_region)

## subset of variants overlapping with exonic regions of protein-coding genes.
d.snv.sub <- d.snvs %>% filter(`Func.refGene` %in% c('exonic', 'exonic;splicing', 'splicing', 'UTR3', 'UTR5', 'UTR5;UTR3')) %>% 
						filter(`ExonicFunc.refGene` %in% c("frameshift insertion", "nonsynonymous SNV", "stopgain", "stoploss")) # 8490 rows
d.cnv.sub <- d.cnvs %>% filter(gene_region %in% c("exonic", "splicing", "UTR3", "UTR5")) # 692 rows
d.sv.sub <- d.svs %>% filter(gene_region %in% c("exonic", "splicing", "UTR3", "UTR5", "UTR5;UTR3")) # 25548 rows

# write.table(d.snv.sub, file.path(data.dir, "case2/case2.snv.final.included.tsv"), quote=F, sep='\t',row.names=F, col.names=T)
# write.table(d.cnv.sub, file.path(data.dir, "case2/case2.cnv.final.included.tsv"), quote=F, sep='\t',row.names=F, col.names=T)
# write.table(d.sv.sub, file.path(data.dir, "case2/case2.sv.final.included.tsv"), quote=F, sep='\t',row.names=F, col.names=T)

## .....

tmp.info <- d.snv.sub$INFO
tmp.depth <- sapply(tmp.info, function(x) stringr::str_extract(x, "DP=\\d+"))
tmp.depth <- stringr::str_replace_all(tmp.depth, "DP=", "") %>% as.numeric()
d.snv.sub$DEPTH = tmp.depth

## compare snvs without read depth filter
genes.with.snvs <- d.snv.sub$`Gene.refGene`
genes.with.snvs <- unique( unlist( sapply(genes.with.snvs, function(x) stringr::str_split(x, ';')[[1]]) ) )

table(genes.with.snvs %in% hscr.gene$Gene)
genes.with.snvs[ genes.with.snvs %in% hscr.gene$Gene ]

## compare snvs with read depth filter
d.snv.sub.dp.filter <- d.snv.sub %>% filter(DEPTH >20)
genes.with.snvs <- d.snv.sub.dp.filter$`Gene.refGene`
genes.with.snvs.list <- unique( unlist( sapply(genes.with.snvs, function(x) stringr::str_split(x, ';')[[1]]) ) )

table(genes.with.snvs.list %in% hscr.gene$Gene)
genes.with.snvs.list[ genes.with.snvs.list %in% hscr.gene$Gene ]

which( grepl("UCHL1", genes.with.snvs) )

d.snv.sub.dp.filter[88:89,]

## ........
## //
## compare cnvs without read depth filter
genes.with.cnvs <- d.cnv.sub$genes
genes.with.cnvs.list <- unique( unlist( sapply(genes.with.cnvs, function(x) stringr::str_split(x, ',')[[1]]) ) )

table(genes.with.cnvs.list %in% hscr.gene$Gene)
genes.with.cnvs.list[ genes.with.cnvs.list %in% hscr.gene$Gene ]

which( grepl( "DCX" , genes.with.cnvs) ) # 459:460
which( grepl( "FGF13" , genes.with.cnvs) ) #   579:582


d.cnv.sub[c(459:460, 579:582),]
write.table(d.cnv, file.path(data.dir, "case2/case2.d.cnv.full.tsv"), quote=F,sep='\t',row.names=F, col.names=T)


##// 
## compare cnvs with read depth filter: DP>20
## no read depth info found in the raw xlsx table.


## ........
## //
## compare svs without read depth filter
d.sv.sub$SU <- sapply(d.sv.sub$info, function(x) stringr::str_split(x, ":")[[1]][2])
d.sv.sub$SU <- as.numeric(d.sv.sub$SU)

genes.with.svs <- d.sv.sub$genes
genes.with.svs.list <- unique( unlist( sapply(genes.with.svs, function(x) stringr::str_split(x, ',')[[1]]) ) ) ## more than 10000 genes were found to have SVs.

write.table(d.sv, file.path(data.dir, "case2/case2.d.sv.full.tsv"), quote=F,sep='\t',row.names=F, col.names=T)

## //
## compare svs with read depth filter
d.sv.sub.dp.filter <- d.sv.sub %>% filter(SU > 20)

genes.with.svs <- d.sv.sub.dp.filter$genes
genes.with.svs.list <- unique( unlist( sapply(genes.with.svs, function(x) stringr::str_split(x, ',')[[1]]) ) ) ## more than 10000 genes were found to have SVs.

table(genes.with.svs.list %in% hscr.gene$Gene)
sort( genes.with.svs.list[ genes.with.svs.list %in% hscr.gene$Gene ] )

which( grepl( "DLX1" , genes.with.svs) ) # 129
which( grepl( "ERBB4" , genes.with.svs) ) #   223:225
which( grepl( "IHH" , genes.with.svs) ) #   242
which( grepl( "TMEFF2" , genes.with.svs) ) #  175
which( grepl( "ZEB2" , genes.with.svs) ) #   67:68


d.sv.sub.dp.filter[c(129, 223:225,242,175,67:68),]
