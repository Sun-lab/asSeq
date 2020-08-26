library(data.table)
library(GenomicRanges)
library(qvalue)
library(readxl)
library(stringr)

specf = "specifications.txt"
getwd()

specs = unlist(read.table(specf, as.is=T))
specs
pref = specs[1]
nsam = specs[2]
queue = specs[3]
days = specs[4]
bmem = as.numeric(specs[5])
seedval = specs[13]
wrk.dir = specs[14]
lib.dir = specs[15]
bas.dir = specs[16]
setwd(wrk.dir)

dyn.dir = sprintf("%s/dynamic", wrk.dir)
if(!file.exists(dyn.dir))dir.create(dyn.dir)


#intro goseq
library("org.Hs.eg.db")

library(GenomeInfoDb)
library(goseq)
library(GenomicFeatures)
head(supportedOrganisms())

get_block = function(x, split="\\.", block=1){
  unlist(strsplit(x, split=split))[block]
}
get_blocks = function(x, split="\\.", blocks=c(1,2,3)){
  unlist(strsplit(x, split=split))[blocks]
}

r38dir = sprintf("%s/Reference", bas.dir)
r38fil = "gencode.v26.GRCh38.genes.gtf"
coni = file(sprintf("%s/c2.cp.reactome.v7.1.symbols.gmt", r38dir))
symb = readLines(coni)
close(coni)
symbl = sapply(symb, strsplit, split="\t")
names(symbl) = 1:length(symbl)
for(i in 1:length(symbl)){
  names(symbl)[i] = symbl[[i]][1]
  symbl[[i]] = symbl[[i]][-(1:2)]
}

symbg = unlist(symbl)
names(symbg)=NULL
syms <- mapIds(org.Hs.eg.db, keys = symbg, column="ENSEMBL", keytype = "SYMBOL", multiVals="first")

varrep = function(i, lst){
  rep(names(lst)[i], each=length(lst[[i]]))
}
symbc = unlist(sapply(1:length(symbl), varrep, lst=symbl))
names(symbc)=NULL
c(length(symbg), length(symbc))



gcats = sprintf("%s/gene2cat_c2.cp.reactome.v7.1.csv", r38dir)
if(!file.exists(gcats)){
  txdb <- makeTxDbFromGFF(sprintf("%s/%s", r38dir, r38fil),format="gtf")
  exons.list.per.gene <- exonsBy(txdb,by="gene")
  # then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
  exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
  exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
  egs = unlist(exonic.gene.sizes); names(egs) = names(exonic.gene.sizes)
  
  coni = file(sprintf("%s/%s", r38dir, r38fil))
  r38inp = readLines(con = coni, n = -1, ok = TRUE, warn = TRUE,
            encoding = "unknown", skipNul = FALSE)
  close(coni)
  length(r38inp)
  tokp = grep("exon", r38inp)
  r38inp = r38inp[tokp]
  length(r38inp)
  r38inp = t(sapply(r38inp, get_blocks, split="\t", blocks=c(1,4,5,7,9)))
  rownames(r38inp) = NULL
  colnames(r38inp) = c("chr", "start", "end", "strand", "info")
  r38inp = data.frame(r38inp,stringsAsFactors = FALSE)
  r38inp[1:4,]
  
  r38inp$id = sapply(r38inp$info, get_blocks, split=";", blocks=1)
  r38inp$nm = sapply(r38inp$info, get_blocks, split=";", blocks=4)
  r38inp[1:4,]
  r38inp$id = gsub("\"", "", gsub("gene_id \"", "", r38inp$id))
  r38inp$nm = gsub("\"", "", gsub(" gene_name \"", "", r38inp$nm))
  r38inp$nm = gsub(" ", "", r38inp$nm)
  
  
  r38inp[1:4,c(-5)]
  #r38inp$len = r38inp$end - r38inp$start
  write.csv(r38inp, sprintf("%s/hg38exons.csv", r38dir), row.names=F, quote=F)
  
  chroms = aggregate(as.character(r38inp[,"chr"]), by=list(r38inp[,"id"], r38inp[,"nm"]), FUN=unique)
  starts = aggregate(as.numeric(as.character(r38inp[,"start"])), by=list(r38inp[,"id"], r38inp[,"nm"]), FUN=min)
  ends =  aggregate(as.numeric(as.character(r38inp[,"end"])), by=list(r38inp[,"id"], r38inp[,"nm"]), FUN=max)
  gninfo = data.frame(chroms)
  gninfo$start =starts[,3]
  gninfo$end=ends[,3]  
  m = match(gninfo[,1], names(egs)); table(names(egs)[m]==gninfo[,1])
  gninfo$lens = egs[m]
  colnames(gninfo)[1:3] = c("id", "nm", "chr")
  gninfo$nm = gsub(" ", "", gninfo$nm)
  gninfo[1:4,]
  write.csv(gninfo, sprintf("%s/hg38genes.csv", r38dir), row.names=F, quote=F)
  gninfo$ids = sapply(gninfo$id, get_block);gninfo[1:2,]
  
  
  m = match(symbg, gninfo$nm);table(is.na(m))
  symbe = gninfo$ids[m]
  gene2cat = data.frame(gene=symbe, category=symbc)
  
  #m = match(syms, gninfo$nm);table(is.na(m))
  #symbe2 = gninfo$ids[m]
  gene2cat2 = data.frame(gene=syms, category=symbc)
  write.csv(gene2cat, sprintf("%s/gene2cat_c2.cp.reactome.v7.1.csv", r38dir), row.names=F)
  write.csv(gene2cat2, sprintf("%s/gene2cat2_c2.cp.reactome.v7.1.csv", r38dir), row.names=F)
}
gninfo = read.csv(sprintf("%s/hg38genes.csv", r38dir), as.is=T)
gninfo[1:2,]
gene2cat2 = read.csv(sprintf("%s/gene2cat2_c2.cp.reactome.v7.1.csv", r38dir), as.is=T)






model = "short"
#model = "long"
for(model in c("short", "long")){
#model = "long"
  cond = "ctcf"
  # --------------------------------------------------------------------
  # read in resutls
  # --------------------------------------------------------------------
  
  diri = "qb_bb2"
  fili = sprintf("%s/%s_%s_%s_final.csv", diri, pref, cond, model)
  resi = fread(fili)
  resi = resi[!is.na(resi$pll.cnd),]
  dim(resi)
  resi$p_cond = resi$pll.cnd
  colnames(resi)[1] = "id"
  resi[1:2,]
  
  
  pi0 = 2*mean(resi$p_cond > 0.5)
  qval = nrow(resi)*pi0*resi$p_cond/rank(resi$p_cond)
  qval = qvalue(resi$p_cond)$qvalue
  table(qval < 0.01)
  table(qval < 0.05)
  table(qval < 0.10)
  table(qval < 0.15)
  table(qval < 0.20)
  table(qval < 0.25)
  
  resi$qval = qval
  
  # --------------------------------------------------------------------
  # read in gene annoation information
  # --------------------------------------------------------------------
  
  gene_file = "gencode.v26.GRCh38.genes_gene_level_anno.txt"
  gene_file = sprintf("%s/Reference/%s", bas.dir, gene_file)
  genes = fread(gene_file)
  dim(genes)
  genes[1:2,]
  
  table(resi$id %in% genes$geneId)
  mat1 = match(resi$id, genes$geneId)
  
  resi = cbind(resi, genes[mat1,])
  dim(resi)
  resi[1:2,]
  
  resi$promoter_start = rep(NA, nrow(resi))
  resi$promoter_end   = rep(NA, nrow(resi))
  
  wn = which(resi$strand == "-")
  wp = which(resi$strand == "+")
  
  resi$promoter_start[wn] = resi$end[wn]
  resi$promoter_end[wn]   = resi$end[wn] + 199
  
  resi$promoter_start[wp] = resi$start[wp] - 199
  resi$promoter_end[wp]   = resi$start[wp]
  
  resi$start.e = resi$start - 1000
  resi$end.e   = resi$end + 1000
  
  resi$promoter_start[wp] = resi$start[wp] - 199
  resi$promoter_end[wp]   = resi$start[wp]
  
  dim(resi)
  resi[1:5,]
  
  gr1 = makeGRangesFromDataFrame(resi, ignore.strand=TRUE, 
                                 seqnames.field="chr",
                                 start.field="promoter_start", 
                                 end.field="promoter_end")
  gr1a = makeGRangesFromDataFrame(resi, ignore.strand=TRUE, 
                                 seqnames.field="chr",
                                 start.field="start.e", 
                                 end.field="end.e")
  
  # --------------------------------------------------------------------
  # read in CTCF binding site information
  # --------------------------------------------------------------------
  
  ff1  = "CTCFBSDB_all_exp_sites_Sept12_2012_hg38_loci.bed.gz"
  ff1 = sprintf("%s/Reference/CTCF/%s", bas.dir, ff1)
  
  resi.bs = read.table(gzfile(ff1), as.is=T)
  dim(resi.bs)
  resi.bs[1:5,]
  
  names(resi.bs) = c("chr", "start", "end")
  table(resi.bs$chr)
  
  lens = resi.bs$end - resi.bs$start + 1
  summary(lens)
  
  pdf(sprintf("%s/Reference/CTCF/CTCFBS_len_hist.pdf", bas.dir), width=6, height=4)
  par(mar=c(5,4,1,1), bty="n")
  hist(log10(lens), xlab="log10(CTCF BS length)", main="", breaks=100)
  abline(v=log10(200))
  dev.off()
  
  table(lens < 400)/length(lens)
  table(lens < 300)/length(lens)
  table(lens < 200)/length(lens)
  
  gr2 = makeGRangesFromDataFrame(resi.bs, ignore.strand=TRUE, 
                                 seqnames.field="chr",
                                 start.field="start", 
                                 end.field="end")
  
  gr3 = makeGRangesFromDataFrame(resi.bs[which(lens < 200),], 
                                 ignore.strand=TRUE, 
                                 seqnames.field="chr",
                                 start.field="start", 
                                 end.field="end")
  
  fun1 <- function(x) sum(width(reduce(x, ignore.strand=T)))
  fun1(gr1)
  
  width2 = fun1(gr2)
  width3 = fun1(gr3)
  width2
  width3
  
  prp2 = width2/(3234.83*10^6)
  prp3 = width3/(3234.83*10^6)
  prp2
  prp3
  
  mtch2 = findOverlaps(gr1, gr2, select="first")
  table(!is.na(mtch2))
  
  mtch3 = findOverlaps(gr1, gr3, select="first")
  table(!is.na(mtch3))
  
  table(qval < 0.25, !is.na(mtch2))
  table(qval < 0.25, !is.na(mtch3))
  
  
  
  qvals = c(0.01, 0.05, 0.10, 0.20, 0.25)
  res2 = matrix(NA, nrow=length(qvals), ncol=5);rownames(res2)=qvals
  colnames(res2) = c("out.qnon", "out.qsig", "in.qnon", "in.qsig", "fisher")
  res3 = res2
  
  for(i in 1:length(qvals)){
    use = sum(qval<qvals[i]);use
    if(use>0){
      tbl = table(qval < qvals[i], !is.na(mtch2))
      res2[i,5] = fisher.test(tbl)$p.value
      res2[i,1:4] = c(tbl)
  
      tbl = table(qval < qvals[i], !is.na(mtch3))
      res3[i,5] = fisher.test(tbl)$p.value
      res3[i,1:4] = c(tbl)
    }
  }
  res2
  res3
  
  out = sprintf("%s/%s_%s_%s_dyn.csv", dyn.dir, pref, cond, model)
  out
  write.csv(res3, out)

  #goseq
  genes = as.integer(qval<qvals[5])
  names(genes) = sapply(resi$id, get_block)
  table(genes)
  #pwf=nullp(genes,"hg19","ensGene")
  m = match(names(genes), gninfo$ids)
  bias.data = gninfo$lens[m]
  pwf=nullp(genes,"hg38",bias.data=bias.data)
  head(pwf)
  
  tokp = rownames(pwf) %in% gene2cat2[,1]; table(tokp)
  path = goseq(pwf, "hg38", gene2cat=gene2cat2)
  
  path[,"over_represented_pvalue"][path[,"over_represented_pvalue"]>1]=1
  path[,"over_represented_pvalue"][path[,"over_represented_pvalue"]<0]=0
  summary(path[,"over_represented_pvalue"])
  qvalue(path[,"over_represented_pvalue"])$pi0
  #table(qvalue(path[,"over_represented_pvalue"])$qvalue<0.25)
  path$qval.over = qvalue(path[,"over_represented_pvalue"])$qvalue
  path$qval.und = qvalue(path[,"under_represented_pvalue"])$qvalue
  over = path[path[,"over_represented_pvalue"]<0.01,]
  under = path[path[,"under_represented_pvalue"]<0.01,]
  out = sprintf("%s/%s_%s_%s_goseq.csv", dyn.dir, pref, cond, model)
  out
  write.csv(rbind(over, under), out, row.names=F)


  
  

  #
  #tp53
  #
  cond = "tp53"
  diri = "qb_bb2"
  fili = sprintf("%s/%s_%s_%s_final.csv", diri, pref, cond, model)
  resi = fread(fili)
  
  dim(resi)
  resi = resi[!is.na(resi$pll.cnd),]
  dim(resi)
  resi[1:2,]
  
  resi$p_cond = resi$pll.cnd
  colnames(resi)[1] = "id"
  
  summary(resi$p_cond)
  pi0 = 2*mean(resi$p_cond > 0.5)
  pi0
  qval = qvalue(resi$p_cond)$qvalue
  table(resi$p_cond < 0.01)
  table(qval < 0.01)
  table(qval < 0.05)
  table(qval < 0.10)
  table(qval < 0.15)
  table(qval < 0.20)
  table(qval < 0.25)
  
  resi$qval = qval
  resi$ens_id = str_extract(resi$id, '(\\S+)(?=\\.\\d+)')
  
  dim(resi)
  resi[1:2,]
  
  # --------------------------------------------------------------------
  # read in TP53 targets
  # --------------------------------------------------------------------
  
  resi_anno = read_excel(sprintf("%s/Reference/TP53/onc2016502x4.xlsx", bas.dir))
  dim(resi_anno)
  resi_anno[1:2,]
  
  table(resi_anno$`ensembl ID` %in% resi$ens_id)
#  resi_anno[!resi_anno$`ensembl ID` %in% resi$ens_id,]
  
  resi_target = resi$ens_id %in% resi_anno$`ensembl ID`
  table(resi_target)
  
#  table(qval < 0.01, resi_target)
#  table(qval < 0.05, resi_target)
#  table(qval < 0.10, resi_target)
#  table(qval < 0.20, resi_target)
  
  qvals = c(0.01, 0.05, 0.10, 0.20, 0.25)
  res2 = matrix(NA, nrow=length(qvals), ncol=5);rownames(res2)=qvals
  colnames(res2) = c("out.qnon", "out.qsig", "in.qnon", "in.qsig", "fisher")
  res3 = res2
  
  for(i in 1:length(qvals)){
    use = sum(qval<qvals[i]);use
    if(use>0){
      tbl = table(qval < qvals[i], resi_target)
      res3[i,5] = fisher.test(tbl)$p.value
      res3[i,1:4] = c(tbl)
    }
  }
  res3
  out = sprintf("%s/%s_%s_%s_dyn.csv", dyn.dir, pref, cond, model)
  out
  write.csv(res3, out)

  #goseq
  genes = as.integer(qval<qvals[5])
  names(genes) = sapply(resi$id, get_block)
  table(genes)
  #pwf=nullp(genes,"hg19","ensGene")
  m = match(names(genes), gninfo$ids)
  bias.data = gninfo$lens[m]
  pwf=nullp(genes,"hg38",bias.data=bias.data)
  head(pwf)
  
  tokp = rownames(pwf) %in% gene2cat2[,1]; table(tokp)
  path = goseq(pwf, "hg38", gene2cat=gene2cat2)
  
  path[,"over_represented_pvalue"][path[,"over_represented_pvalue"]>1]=1
  path[,"over_represented_pvalue"][path[,"over_represented_pvalue"]<0]=0
  summary(path[,"over_represented_pvalue"])
  qvalue(path[,"over_represented_pvalue"])$pi0
  #table(qvalue(path[,"over_represented_pvalue"])$qvalue<0.25)
  path$qval.over = qvalue(path[,"over_represented_pvalue"])$qvalue
  path$qval.und = qvalue(path[,"under_represented_pvalue"])$qvalue
  over = path[path[,"over_represented_pvalue"]<0.01,]
  under = path[path[,"under_represented_pvalue"]<0.01,]
  out = sprintf("%s/%s_%s_%s_goseq.csv", dyn.dir, pref, cond, model)
  out
  write.csv(rbind(over, under), out, row.names=F)


  #age
  cond = "age"
  diri = "qb_bb2"
  fili = sprintf("%s/%s_%s_%s_final.csv", diri, pref, cond, model)
  resi = fread(fili)
  
  dim(resi)
  resi = resi[!is.na(resi$pll.cnd),]
  dim(resi)
  resi[1:2,]
  
  resi$p_cond = resi$pll.cnd
  colnames(resi)[1] = "id"
  
  summary(resi$p_cond)
  pi0 = 2*mean(resi$p_cond > 0.5)
  pi0
  qval = qvalue(resi$p_cond)$qvalue
  table(resi$p_cond < 0.01)
  table(qval < 0.01)
  table(qval < 0.05)
  table(qval < 0.10)
  table(qval < 0.15)
  table(qval < 0.20)
  table(qval < 0.25)
}





q(save = "no")


