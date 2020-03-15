args = commandArgs(trailingOnly = TRUE)
#args = c("21", "1", "194", "100", "5e5", "long")
#args = c("22", "1", "194", "100", "5e5", "long")
args
chri = as.numeric(args[1])
blocki = as.numeric(args[2])
nsub = as.numeric(args[3])
numpoints = as.numeric(args[4])
cis_window = as.numeric(args[5])
model = args[6]
suff0 = sprintf("%s_%s", chri, blocki)

#numpoints = 50
library(MatrixEQTL)

useModel = modelLINEAR; 
set.seed(12345)
#nsub = 280
#cisDist = 5e5
#model = "long"
pref = "Whole_Blood"


#root.dir = "/pine/scr/z/h/zhabotyn/R01"
#root.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/_GTEx"
#work.dir = sprintf("%s/2019_10_10", root.dir)
#setwd(work.dir)
library(MatrixEQTL)
useModel = modelLINEAR; 
source("helpers.R")

out.dir = sprintf("oneperm_%s_%s_%s_%s", pref, nsub, cis_window, model)
perm.dir = sprintf("boot_%s_%s_%s_%s_%s", pref, nsub, cis_window, model, numpoints)

#boot.dir = sprintf("bootSub%s", nsub)
int.dir = sprintf("%s_%s_%s", pref, nsub, cis_window)


  
  genotype_file_name = sprintf("%s/genotypes_%s.dat", int.dir, suff0)
  covariates_file_name = sprintf("%s/Xmat_%s.csv", int.dir, model)
  covar = read.csv(covariates_file_name,header=F,as.is=T)
#  SNP_file_name = sprintf("%s/SNP_%s.txt", int.dir, suff0)
#  covariates_file_name = sprintf("%s/Covariates.txt", boot.dir)  
#  covariates_file_name = sprintf("%s/Covariates.txt", boot.dir)  
#  covar = read.table(covariates_file_name,header=T,as.is=T)

#  nms = covar[,1]
#  covar = as.matrix(covar[,-1])
  covar = as.matrix(covar)
#  rownames(covar) = nms
  cvrt = SlicedData$new()
  cvrt = cvrt$CreateFromMatrix(t(covar))

#  gen.in = read.table(SNP_file_name, header=T)
  gen.in = read.table(genotype_file_name, header=T)
  gen.in[gen.in==3] = 1
  gen.in[gen.in==4] = 2
  
#  snpspos_file_name = sprintf("%s/snpspos_%s.dat", int.dir, suff0)
  snpspos_file_name = sprintf("%s/genotypei_%s.dat", int.dir, suff0)
  snpspos = read.table(file=snpspos_file_name, header=T)
  #nms = gen.in[,1]
  #gen.in = as.matrix(gen.in[,-1])
  gen.in = as.matrix(gen.in)
  rownames(gen.in) = snpspos[,1]
  #rownames(gen.in) = nms
  snps = SlicedData$new();
  snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
  snps = snps$CreateFromMatrix(gen.in)

  expression_file_name = sprintf("%s/GE_norm_%s_%s.dat", int.dir, model, suff0)
  exprj = read.table(expression_file_name, header=F)
#  exprj = read.csv(expression_file_name, header=F)

  output_file_name2 = sprintf("%s/output_eigenMT_%s.txt", out.dir, suff0)

  gen.sub = read.csv(sprintf("%s/min_snp_vals_%s.csv", out.dir, suff0))
  eigenMT = read.csv(sprintf("%s/upd_eigenMT_%s.csv", out.dir, suff0))

  genepos_file_name = sprintf("%s/genepos_%s.dat", int.dir, suff0)
  genepos = read.table(file=genepos_file_name, header=T)

  pvOutputThreshold = 1;
  errorCovariance = numeric();

#
#
#produce boots
#
#


#first run trial to see how much lower boots are wrt permu p-vals
#then run all 1000 permutations
#produce
target.perm.ps = 10^-seq(-log10(0.001), -log10(0.5), length.out=numpoints)
target.perm.pi = target.perm.ps[1]
eigenMTp = eigenMT
j = 0
tmpres = matrix(NA, nrow=10, ncol=3)
repeat{
  j = j + 1
  redboot = t(sapply(1:length(target.perm.ps), get_reduced_boot, target.perm.ps=target.perm.ps, i=1, 
  mQTL.fit=eigenMTp, expr.mat = exprj, min.SNP=gen.sub, covars=covar, nsam=nsub, use.norm=F))
  rownames(redboot) = sprintf("%s_%s", rownames(exprj)[1], 1:length(target.perm.ps))

  geneb = SlicedData$new();
  geneb = geneb$CreateFromMatrix(redboot)
  output_file_name = sprintf("%s/output_boot_%s.txt", perm.dir, suff0)
  bootpos = genepos[rep(1, nrow(redboot)),]
  bootpos[,1] = rownames(redboot)
  tstarts = proc.time()
  #for(i in 1:1e1){
  meb = Matrix_eQTL_main(
        snps = snps,
        gene = geneb,
        cvrt = cvrt,
        pvOutputThreshold = 1e-200,
        output_file_name = sprintf("%s_tmp", output_file_name),
        output_file_name.cis = output_file_name,
        pvOutputThreshold.cis = 1e-200,
        useModel = useModel, 
        errorCovariance = errorCovariance,
        snpspos = snpspos,
        genepos = bootpos, 
        cisDist = 1e9,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = TRUE,
        noFDRsaveMemory = TRUE);
  #}
  tends = proc.time()
  tends[3]-tstarts[3]
  #pvalb = aggregate(meb$cis$eqtls$pvalue, by=list(meb$cis$eqtls$gene), FUN=min)
  pvalb = meb$cis$min.pv.gene
  pvalb = data.frame(names(pvalb), pvalb, stringsAsFactors=FALSE)  
  gns = sort(unique(pvalb[,1]))
  #ords = data.frame(t(sapply(gns, minord, genes=meb$cis$eqtls$gene, pvals=meb$cis$eqtl$pvalue)))
  #ords[,2] = as.numeric(as.character(ords[,2]))
  #ords$betas=meb$cis$eqtl$beta[ords[,2]]
  #ords$tstat=meb$cis$eqtls$statistic[ords[,2]]
  #ords$pvals=meb$cis$eqtls$pvalue[ords[,2]]
  #ords$snps=meb$cis$eqtls$snps[ords[,2]]
  #ords
  ords = pvalb

nperm = 1e2
permind = 1:nperm
permm = matrix(NA, nrow=nperm, ncol=nrow(pvalb))
#do 1k perm
starts = proc.time()
for(i in permind){
  permi = sample(1:nsub)
  exprp = redboot[,permi]
  covap = covar[permi,]
  output_file_name = sprintf("%s/output_boot_%s.txt", perm.dir, suff0)

  cvrtp = SlicedData$new()
  cvrtp = cvrtp$CreateFromMatrix(t(covap))
    
  genep = SlicedData$new();
  genep = geneb$CreateFromMatrix(exprp)
  
  mep = Matrix_eQTL_main(
          snps = snps,
          gene = genep,
          cvrt = cvrtp,
          pvOutputThreshold = 1e-200,
          output_file_name = sprintf("%s_tmp", output_file_name),#permout_file_name),
          output_file_name.cis = output_file_name,#permout_file_name,
          pvOutputThreshold.cis = 1e-200,
          useModel = useModel, 
          errorCovariance = errorCovariance,
          snpspos = snpspos,
          genepos = bootpos, 
          cisDist = 1e9,
          verbose = TRUE,
          pvalue.hist = TRUE,
          min.pv.by.genesnp = TRUE,
          noFDRsaveMemory = TRUE);
  #pvalp = aggregate(mep$cis$eqtls$pvalue, by=list(mep$cis$eqtls$gene), FUN=min)
  pvalp = mep$cis$min.pv.gene
  pvalp = data.frame(names(pvalp), pvalp, stringsAsFactors=FALSE)  
  permm[i,] = pvalp[,2]
  message(i, " out of ", nperm)
}
ends = proc.time()
tim0 = (ends[3]-starts[3])
message("time per iter: ", tim0/nperm)#half minute per iteration per 1 gene
colnames(permm) = pvalb[,1]
#50 points - 1 second (for chr 22) //
y = rep(NA, ncol(permm))
for(i in 1:ncol(permm)){
  y[i] = sum(permm[1:nperm,i]< pvalb[i,2])
}
message("attempt ", j, ": ", mean(y<(nperm*0.001)), " and ", mean(y>(nperm*.3)))
tmpres[j, 1] = sum(y<(nperm*0.001))
tmpres[j, 2] = sum(y>(nperm*.3))
tmpres[j, 3] = tim0

  if(mean(y<(nperm*0.001))>.40 & j<5){
    #consider the case with all 0's
    if(eigenMTp$TESTS<=eigenMT$TESTS){
      eigenMTp$TESTS = eigenMTp$TESTS/2
    }else{
      eigenMTp$TESTS = eigenMTp$TESTS - abs(eigenMT$TESTS-eigenMTp$TESTS)/2
    }
  }else{
    if(mean(y>(nperm*.3))>.3 & j<5){
      if(eigenMTp$TESTS>=eigenMT$TESTS){
        eigenMTp$TESTS = eigenMTp$TESTS*2      
      }else{
        eigenMTp$TESTS = eigenMTp$TESTS + abs(eigenMT$TESTS-eigenMTp$TESTS)/2    
      }
    }else{
      break
    }
  }
}
tmpres
#
#the rest of iterations
#
nperm0 = nperm
nperm = 1e3
permind = (nperm0+1):nperm
permm0 = permm
permm = matrix(NA, nrow=nperm, ncol=nrow(pvalb))
permm[1:nperm0,] = permm0


starts = proc.time()
for(i in permind){
  permi = sample(1:nsub)
  exprp = redboot[,permi]
  covap = covar[permi,]
  output_file_name = sprintf("%s/output_boot_%s.txt", perm.dir, suff0)

  cvrtp = SlicedData$new()
  cvrtp = cvrtp$CreateFromMatrix(t(covap))
    
  genep = SlicedData$new();
  genep = geneb$CreateFromMatrix(exprp)
  
  mep = Matrix_eQTL_main(
          snps = snps,
          gene = genep,
          cvrt = cvrtp,
          pvOutputThreshold = 1e-200,
          output_file_name = sprintf("%s_tmp", output_file_name),#permout_file_name),
          output_file_name.cis = output_file_name,#permout_file_name,
          pvOutputThreshold.cis = 1e-200,
          useModel = useModel, 
          errorCovariance = errorCovariance,
          snpspos = snpspos,
          genepos = bootpos, 
          cisDist = 1e9,
          verbose = TRUE,
          pvalue.hist = TRUE,
          min.pv.by.genesnp = TRUE,
          noFDRsaveMemory = TRUE);
  #pvalp = aggregate(mep$cis$eqtls$pvalue, by=list(mep$cis$eqtls$gene), FUN=min)
  pvalp = mep$cis$min.pv.gene
  pvalp = data.frame(names(pvalp), pvalp, stringsAsFactors=FALSE)  
  permm[i,] = pvalp[,2]
  message(i, " out of ", nperm)
}
ends = proc.time()
tim1 = (ends[3]-starts[3])
message("time per iter: ", (tim0+tim1)/nperm)#half minute per iteration per 1 gene



y = rep(NA, ncol(permm))
for(i in 1:ncol(permm)){
  y[i] = sum(permm[1:nperm,i]< pvalb[i,2])
}
table(y)
message("attempt ", j, ": ", mean(y<(nperm*0.001)), " and ", mean(y>nperm*.3))


#kp0 = (y/nperm)>0.002 & (y/nperm)<.25;table(kp0)
kp3 = (y/nperm)>=0     & (y/nperm)<=0.3
kp3a = (y/nperm)>0     & (y/nperm)<=0.3
    
y1 = -log10(y/nperm)
x1 = -log10(pvalb[,2])
glmi3 = glm(cbind(y[kp3],nperm-y[kp3])~x1[kp3], family="binomial")
#summary(glmi3)
lmi3 = lm(y1[kp3a]~x1[kp3a])
#summary(lmi3)


eigenMT$TESTSupd = eigenMTp$TESTS
xval = -log10(eigenMT$p.value)
pred.perm = logiti(glmi3$coef[1]+glmi3$coef[2]*xval)
pred.perm
eigenMT$pred.permGLM = pred.perm
pred.perm = 10^-(lmi3$coef[1] +lmi3$coef[2]*xval)
pred.perm[pred.perm>1] = 1
pred.perm
eigenMT$pred.permLM = pred.perm

eigenMT$LM.i = lmi3$coef[1]
eigenMT$LM.s = lmi3$coef[2]
eigenMT$GLM.i = glmi3$coef[1]
eigenMT$GLM.s = glmi3$coef[2]
eigenMT$numpts = sum(kp3a)

#write.csv(eigenMT, sprintf("%s/upd_eigenMT_%s.csv", work.dir, suff0), quote=F, row.names=F)
write.csv(eigenMT, sprintf("%s/upd_eigenMT_%s.csv", perm.dir, suff0), quote=F, row.names=F)

filout = sprintf("%s/boot_pval_%s.csv", perm.dir, suff0)
filout
write.csv(permm, filout, row.names=F, quote=F)
#perm.dir
permm2 = rbind(permm, pvalb[,2])
#ords$pvals,
#ords$betas,
#ords$tstat)
colnames(permm2) = ords[,1]

write.csv(permm2, filout, row.names=F, quote=F)

permp = rep(NA, ncol(permm))
for(i in 1:ncol(permm)){
  permp[i] = mean(permm2[nperm+1,i]>=permm[1:nperm,i])
}
permp
filout = sprintf("%s/short_boot_pval_%s.csv", perm.dir, suff0)
filout
write.csv(cbind(pvalb, permp), filout, row.names=F, quote=F)

filout = sprintf("%s/time_%s.csv", perm.dir, suff0)
filout
write.csv(tim0+tim1, filout, row.names=F, quote=F)

file.remove(sprintf("%s_tmp", output_file_name))
file.remove(sprintf("%s", output_file_name))

q("no")
