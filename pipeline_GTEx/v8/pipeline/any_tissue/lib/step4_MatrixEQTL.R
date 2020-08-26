#sbatch -p general -N 1 -t 01-00:00:00 -o batch.out --mem=32g --wragenepos_p="R CMD BATCH step3_MatrixEQTL.R step3_MatrixEQTL.Rout"
args = commandArgs(trailingOnly = TRUE)
#need to fix chromosomes 7 and 8
args
chri = as.numeric(args[1])
nsub = as.numeric(args[2])
seedval = as.numeric(args[3])
cis_window = as.numeric(args[4])
model = args[5]
if(length(args)>5){
  paral = as.numeric(args[6])
}else{
  paral = 1e6
}

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
#eigenMTdir = "/nas/longleaf/home/zhabotyn/pythlab/eigenMT"
eigenMTdir = spects[9]
setwd(wrk.dir)

library(MatrixEQTL)
library(Matrix)
useModel = modelLINEAR; 
source(sprintf("%s/helpers.R", lib.dir))


numpoints = 100
maf = 0.05
        
set.seed(seedval)


#norm = F
#norm = T
#data.dir = "../data"
#cnt.dir = sprintf("%s/cnt", data.dir)
#geno.dir = "../datagen" 

routdir = sprintf("%s/rout_%s", wrk.dir, pref)
boutdir = sprintf("%s/bout_%s", wrk.dir, pref)
if(!file.exists(routdir))dir.create(routdir)
if(!file.exists(boutdir))dir.create(boutdir)

int.dir = sprintf("%s_%s_%s", pref, nsub, cis_window)

#root.dir = "/pine/scr/z/h/zhabotyn/R01"
#gtex.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8"
#cnt.dir = sprintf("%s/%s_prepr", gtex.dir, pref)
cnt.dir = sprintf("%s_prepr", pref);cnt.dir;file.exists(cnt.dir)
out.dir = sprintf("oneperm_%s_%s_%s_%s", pref, nsub, cis_window, model)
perm.dir = sprintf("boot_%s_%s_%s_%s_%s", pref, nsub, cis_window, model, numpoints)
if(!file.exists(out.dir))dir.create(out.dir)
if(!file.exists(perm.dir))dir.create(perm.dir)


genepos_file_name = sprintf("%s/geneInfo_prepr_%s.txt", cnt.dir, model)
geneInfo = read.table(genepos_file_name, 
                      header = T, as.is = T)
genepos = geneInfo[geneInfo$chr==sprintf("chr%s", chri),1:4]
genepos[,2] = gsub("chr", "", genepos[,2])
for(coli in 3:4)genepos[,coli] = as.numeric(genepos[,coli])


covariates_file_name = sprintf("%s/Xmat_%s.csv", int.dir, model) 
covar =  read.csv(covariates_file_name, as.is=T, header=F)
covar = as.matrix(covar)

converge = 1e-4
vari = apply(covar,2,var)

updvar = which(vari<converge)
for(i in updvar){
  if(length(vari[-updvar]>0)>0){
    correct = sqrt(median(vari[-updvar]))/sqrt(vari[i])
  }else{
    correct = 1/sqrt(vari[i])
  }    
  xm = mean(covar[,i])    
  covar[,i] = xm+(covar[,i]-xm)*correct
}
vari
apply(covar,2,var)

  blocki = 1
  countjobs = 0

for(blocki in 1:nrow(genepos)){
#  for(blocki in 1:10){
  suff0 = sprintf("%s_%s", chri, blocki)
  timout = sprintf("%s/time_%s.csv", perm.dir, suff0)
  if(!file.exists(timout)){
    results = tryCatch({
    #suff0 = sprintf("%s_%s", chri, blocki)
    output_file_name = sprintf("%s/output_norm_%s.txt", int.dir, suff0)
    output_file_name2 = sprintf("%s/output_eigenMT_%s.txt", out.dir, suff0)
    expression_file_name = sprintf("%s/GE_norm_%s.dat", int.dir, suff0)
    output_file_name_min = sprintf("%s/output_norm_min_%s.txt", perm.dir, suff0)


    genotype_file_name = sprintf("%s/genotypes_%s.dat", int.dir, suff0)
    cvrt = SlicedData$new()
    cvrt = cvrt$CreateFromMatrix(t(covar))
  
    g.ini = read.table(genotype_file_name, header=T)
    g.ini[g.ini==3] = 1
    g.ini[g.ini==4] = 2
    snpspos_file_name = sprintf("%s/genotypei_%s.dat", int.dir, suff0)
    snpspos = read.table(snpspos_file_name, header=T, as.is=T)
    for(coli in 3:3)snpspos[,coli] = as.numeric(snpspos[,coli])
    rownames(g.ini) = snpspos[,1]

    kp = rowMeans(g.ini)/2
    
    converge=5e-5
    varZ = apply(g.ini, 1, var)
    wVar = (varZ >= converge)
    kp = wVar #& ((a0&a1)|(a2&a1)|(a0&a2))
    table(kp)

    SNP_file_name = sprintf("%s/SNP_%s.txt", int.dir, suff0)

    write.table(g.ini, SNP_file_name, row.names=T, col.names=T, quote=F, sep="\t")
  
    expression_file_name = sprintf("%s/GE_norm_%s_%s.dat", int.dir, model, suff0)
    exprj = read.table(expression_file_name)
    
    pvOutputThreshold = 1;
    errorCovariance = numeric();
      
    snps = SlicedData$new();
    snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
    snps = snps$CreateFromMatrix(as.matrix(g.ini))
      
    genepos_file_name = sprintf("%s/genepos_%s.dat", int.dir, suff0)
    colnames(snpspos) = c("snpid", "chr", "pos")
    colnames(genepos) = c("geneid", "chr", "left", "right")
    write.table(genepos[blocki,], file=genepos_file_name, row.names=F, col.names=T, quote=F, sep="\t")


    rownames(exprj) = genepos$geneid[blocki]
    gene = SlicedData$new();
    gene = gene$CreateFromMatrix(as.matrix(exprj))

    me = Matrix_eQTL_main(
          snps = snps,
          gene = gene,
          cvrt = cvrt,
          pvOutputThreshold = 1e-200,
          output_file_name = sprintf("%s_tmp", output_file_name),
          output_file_name.cis = output_file_name,
          pvOutputThreshold.cis = pvOutputThreshold,
          useModel = useModel, 
          errorCovariance = errorCovariance,
          snpspos = snpspos,
          genepos = genepos[blocki,], 
          cisDist = 1e9,
          verbose = TRUE,
          pvalue.hist = TRUE,
          min.pv.by.genesnp = FALSE,
          noFDRsaveMemory = FALSE);
  
    file.remove(sprintf("%s_tmp", output_file_name))
    
    #run eigenMT
    cmdi = sprintf("python %s/eigenMT.py --CHROM %s --QTL %s --GEN %s --GENPOS %s --PHEPOS %s --OUT %s", 
                    eigenMTdir, chri, output_file_name, SNP_file_name, snpspos_file_name, 
                     genepos_file_name, output_file_name2)
    system(cmdi)
    
    eigenMT = read.table(output_file_name2, header=T, as.is=T)
    m = match(eigenMT$SNP, snpspos$snpid)
    eigenMT$chr = chri
    eigenMT$snppos = snpspos$pos[m]
    eigenMT$genestart = genepos$left[blocki]
    eigenMT$geneend = genepos$right[blocki]
    #get minimum p-values and respective snps
    genes=as.character(me$cis$eqtls$gene)
    ords = data.frame(t(sapply(sort(unique(genes)), minord, genes=me$cis$eqtls$gene, pvals=me$cis$eqtl$pvalue)))
    ords[,2] = as.numeric(as.character(ords[,2]))
    pvals = aggregate(me$cis$eqtls$pvalue, by=list(me$cis$eqtls$gene), FUN=min)
    table(pvals[,2]==me$cis$eqtl$pvalue[ords[,2]])
    ords$betas=me$cis$eqtl$beta[ords[,2]]
    ords$tstat=me$cis$eqtls$statistic[ords[,2]]
    ords$pvals=me$cis$eqtls$pvalue[ords[,2]]
    ords$snps=me$cis$eqtls$snps[ords[,2]]
    ords
    ntest = aggregate(rep(1, length(me$cis$eqtls$pvalue)), by=list(me$cis$eqtls$gene), FUN=sum)
    m = match(ords[,1], ntest[,1])
    table(ntest[m,1]==ords[,1])
    ords$ntest = ntest[m,2]
    
    nmedp = aggregate(me$cis$eqtls$pvalue, by=list(me$cis$eqtls$gene), FUN=median)
    m = match(ords[,1], nmedp[,1])
    table(nmedp[m,1]==ords[,1])
    ords$nmedp = nmedp[m,2]
    
    m = match(ords[,1], eigenMT$gene)
    m
    table(ords[,1] == eigenMT$gene[m])
    ords$TESTS = eigenMT$TESTS[m]
    eigenMT$ntest[m] = ords$ntest
    eigenMT = eigenMT[m,]
    
    #need to refit minimums with linear model to get other covariates
    m = match(eigenMT$SNP, rownames(g.ini))
    table(eigenMT$SNP==rownames(g.ini)[m])
    gen.sub = matrix(g.ini[m,],nrow=length(m));rownames(gen.sub) = rownames(g.ini)[m]
    colnames(gen.sub) = colnames(g.ini)
    write.csv(gen.sub, sprintf("%s/min_snp_vals_%s.csv", out.dir, suff0), quote=F, row.names=F)
    table(rownames(gen.sub)==eigenMT$SNP)
    write.csv(eigenMT, sprintf("%s/upd_eigenMT_%s.csv", out.dir, suff0), quote=F, row.names=F)
    
    
    filout = sprintf("%s/short_pval_%s.csv", out.dir, suff0)
    write.table(ords[,-c(2)], filout, sep=",", row.names=F, col.names=T)
    
    
    
    #write intermediate objects
    write.csv(exprj, sprintf("%s/expr_%s.csv", out.dir, suff0), quote=F)
    write.csv(gen.sub, sprintf("%s/msnp_%s.csv", out.dir, suff0), quote=F)
  
    rinpdir = lib.dir
    rinp = sprintf("%s/step4_runboot.R", rinpdir)
    rout = sprintf("%s/step4_runboot_%s_%s_%s_%s_%s_%s_%s.out",
                    routdir, chri, blocki, nsub, numpoints, cis_window, model, seedval) 
    qout = sprintf("%s/step4_runboot_%s_%s_%s_%s_%s_%s_%s.Rout",
                    boutdir, chri, blocki, nsub, numpoints, cis_window, model, seedval) 
    com = sprintf("R CMD BATCH '--args %s %s %s %s %s %s' %s %s",
                    chri, blocki, numpoints, cis_window, model, seedval, rinp, rout)
    com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                              queue, days, qout, mem, com)        
    message(com2)
    if(blocki%%paral==0){
      system(com)
    }else{
      system(com2)
    }
              countjobs = countjobs + 1
    0
    }, error = function(e){
    1
    })
    if(results==1)message(blocki)
  }
}
message(countjobs)

q("no")


