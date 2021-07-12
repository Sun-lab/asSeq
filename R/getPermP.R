require(deSeq2)
require(MatrixEQTL)

get_blocks = function(x, split="\\|", blocks=c(1)){
  paste(unlist(strsplit(x, split=split))[blocks], collapse="-")
}
get_blockn = function(x, split="_", blocks=2:3){
  unlist(strsplit(x, split=split))[blocks]
}

splitting = function(x, split="\\|", block=5, convert=T){
  out = unlist(strsplit(x, split=split))[block]
  if(convert)out = as.numeric(out)
  out
}


minord = function(genei, pvals, genes){
  kp = which(genes==genei)
  c(as.character(genei), kp[order(pvals[kp])[1]])
}

classify = function(i, frac1, fracN, kp){
  out = c(mean(abs(frac1[i,kp[i,]]-.5)>.4, na.rm=T), mean(fracN[i,kp[i,]], na.rm=T), sum(kp[i,], na.rm=T))
}



normscore = function(vec) {
    len  = length(na.omit(vec))+1
    rank = rank(na.omit(vec))
    ties = (rank - floor(rank)) > 0
    new.vec = vec[!is.na(vec)] 
    new.vec[!ties]=qnorm(rank[!ties]/len)
    new.vec[ties] =0.5*(qnorm((rank[ties]+0.5)/len)+qnorm((rank[ties]-0.5)/len))
    vec[!is.na(vec)] = new.vec
    vec
}

logit = function(x){log(x)-log(1-x)}
logiti = function(x)1/(1+exp(-x))
updminx = function(y, a, b){(-log(-1+1/y)-a)/b}

replaceOutliers <- function(object, trim=.2, cooksCutoff, minReplicates=1, whichSamples) {
  if (is.null(attr(object,"modelMatrix")) | !("cooks" %in% assayNames(object))) {
    stop("first run DESeq, nbinomWaldTest, or nbinomLRT to identify outliers")
  }
  if (minReplicates < 1) {
    stop("at least 3 replicates are necessary in order to indentify a sample as a count outlier")
  }
  stopifnot(is.numeric(minReplicates) & length(minReplicates) == 1)
  p <- ncol(attr(object,"modelMatrix"))
  m <- ncol(object)
  if (m <= p) {
    assays(object)[["originalCounts"]] <- counts(object)
    return(object)
  }
  if (missing(cooksCutoff)) {
    cooksCutoff <- qf(.99, p, m - p)
  }
  idx <- which(assays(object)[["cooks"]] > cooksCutoff)
  mcols(object)$replace <- apply(assays(object)[["cooks"]], 1, function(row) any(row > cooksCutoff))
  mcols(mcols(object),use.names=TRUE)["replace",] <- DataFrame(type="intermediate",description="had counts replaced")
  trimBaseMean <- apply(counts(object,normalized=TRUE),1,mean,trim=trim)
  # build a matrix of counts based on the trimmed mean and the size factors
  replacementCounts <- if (!is.null(normalizationFactors(object))) {
    as.integer(matrix(rep(trimBaseMean,ncol(object)),ncol=ncol(object)) * 
                 normalizationFactors(object))
  } else {
    as.integer(outer(trimBaseMean, sizeFactors(object), "*"))
  }
  # replace only those values which fall above the cutoff on Cook's distance
  newCounts <- counts(object)
  newCounts[idx] <- replacementCounts[idx]
  
  if (missing(whichSamples)) {
    whichSamples <- nOrMoreInCell(attr(object,"modelMatrix"), n = minReplicates)
  }
  stopifnot(is.logical(whichSamples))
  object$replaceable <- whichSamples
  mcols(colData(object),use.names=TRUE)["replaceable",] <- DataFrame(type="intermediate",
                                                                     description="outliers can be replaced")
  assays(object)[["originalCounts"]] <- counts(object)
  if (sum(whichSamples) == 0) {
    return(object)
  }
  counts(object)[,whichSamples] <- newCounts[,whichSamples,drop=FALSE]
  object
}

nOrMoreInCell <- function(modelMatrix, n) {
  numEqual <- sapply(seq_len(nrow(modelMatrix)), function(i) {
    modelMatrixDiff <- t(t(modelMatrix) - modelMatrix[i,])
    sum(apply(modelMatrixDiff, 1, function(row) all(row == 0)))
  })
  numEqual >= n
}



#setClass("permEst", slots=list(snpM="numeric", genem="numeric", cvrtm="numeric", 
#pvOutputThreshold="numeric",  pvOutputThreshold.csv="numeric", )



#permEst
#list(snpM, geneM, cvrtM, 
#snpspos, genepos, outpf, 
#pvOutputThreshold=1e-1000, pvOutputThreshold.csv=1, cisDist=1e9,
#effNumGuess=nrow(snpM), 
#verbose=FALSE, pvalue.hist=FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory=FALSE,
#updNtests=NA)
getPermP = function(permEst, num.grid = 100, slice=2000, n.perm=1000, ini.perm=100, nsub=NA, write.res=F, max.adj=10){ 
  #num.grid = 100; slice=2000; n.perm=1000; ini.perm=100; nsub = NA

  require(MatrixEQTL)
  useModel=modelLINEAR
  errorCovariance = numeric();
  if(is.na(nsub))nsub=ncol(permEst$geneM);nsub
  snps = SlicedData$new();
  snps$fileSliceSize = slice
  snps = snps$CreateFromMatrix(permEst$snpM)
      
  gene = SlicedData$new();
  gene = gene$CreateFromMatrix(permEst$geneM)

  cvrt = SlicedData$new()
  cvrt = cvrt$CreateFromMatrix(t(permEst$cvrtM))

  bootfil = sprintf("boot_%s", permEst$outpf)
  if(!is.na(permEst$outdir)){
    if(!file.exists(permEst$outdir))dir.create(permEst$outdir)
    outpf = sprintf("%s/%s", permEst$outdir, permEst$outpf)
    bootfil = sprintf("%s/%s", permEst$outdir, bootfil)
  }
  tmpfil = sprintf("%s_tmp", permEst$outpf)
  tmpbootfil = sprintf("%s_tmp", bootfil)

  me = Matrix_eQTL_main(
          snps = snps,
          gene = gene,
          cvrt = cvrt,
          pvOutputThreshold = permEst$pvOutputThreshold,
          output_file_name = tmpfil,
          output_file_name.cis = permEst$outpf,
          pvOutputThreshold.cis = permEst$pvOutputThreshold.csv,
          useModel = useModel, 
          errorCovariance = errorCovariance,
          snpspos = permEst$snpspos,
          genepos = permEst$genepos, 
          cisDist = permEst$cisDist,
          verbose = permEst$verbose,
          pvalue.hist = permEst$pvalue.hist,
          min.pv.by.genesnp = permEst$min.pv.by.genesnp,
          noFDRsaveMemory = permEst$noFDRsaveMemory);
  if(file.exists(tmpfil))file.remove(tmpfil)
  names(me)
  names(me$cis)


  #
  #
  #produce boots
  #
  #

  #find minimum SNP from me
  #get it's values
  o = order(me$cis$eqtls$pvalue)
  m = match(me$cis$eqtls$snps[o[1]], rownames(permEst$snpM));m
  gen.sub = permEst$snpM[m,,drop=F]
  eigenMTp = me$cis$eqtls[o[1],]
  colnames(eigenMTp) = c("SNP", "gene", "t.stat", "p.value", "FDR", "beta")
  eigenMTp = eigenMTp[,c("SNP", "gene", "beta", "t.stat", "p.value", "FDR")]
  eigenMTp$BF = NA
  eigenMTp$TESTS = permEst$effNumGuess
  eigenMTp$BF = eigenMTp$TESTS*eigenMTp$p.value
  eigenMTp$chr = permEst$genepos[1,2] 
  eigenMTp$snppos = permEst$snpspos[m,3] 
  eigenMTp$genestart = permEst$genepos[1,3] 
  eigenMTp$geneend = permEst$genepos[1,4] 
  eigenMTp$ntest = nrow(permEst$snpM)  
  eigenMT = eigenMTp
 
  #first run trial to see how much lower boots are wrt permu p-vals
  #then run all 1000 permutations
  #produce
  target.perm.ps = 10^-seq(-log10(0.001), -log10(0.5), length.out=num.grid)
  target.perm.pi = target.perm.ps[1]
  eigenMTp = eigenMT
  j = 0
  tmpres = matrix(NA, nrow=100, ncol=4)

  #another update: don't discard iterations that were within range, use them
  repeat{
    j = j + 1                                                                                        
    redboot = t(sapply(1:length(target.perm.ps), get_reduced_boot, target.perm.ps=target.perm.ps, 
    mQTL.fit=eigenMTp, expr.mat = permEst$geneM, min.SNP=gen.sub, covars=permEst$cvrtM, nsam=nsub, use.norm=F))
    if(j>1){
      redboot = rbind(redboot, kpboot)
    }
    rownames(redboot) = sprintf("%s_%s", rownames(permEst$geneM)[1], 1:nrow(redboot))
  
    geneb = SlicedData$new();
    geneb = geneb$CreateFromMatrix(redboot)
    bootpos = permEst$genepos[rep(1, nrow(redboot)),]
    bootpos[,1] = rownames(redboot)
    tstarts = proc.time()
    meb = Matrix_eQTL_main(
          snps = snps,
          gene = geneb,
          cvrt = cvrt,
          pvOutputThreshold = permEst$pvOutputThreshold,
          output_file_name = tmpbootfil,
          output_file_name.cis = bootfil,
          pvOutputThreshold.cis = permEst$pvOutputThreshold,
          useModel = useModel, 
          errorCovariance = errorCovariance,
          snpspos = permEst$snpspos,
          genepos = bootpos, 
          cisDist = permEst$cisDist,
          verbose = permEst$verbose,
          pvalue.hist = permEst$verbose,
          min.pv.by.genesnp = TRUE,
          noFDRsaveMemory = TRUE);
    #}
    tends = proc.time()
    tends[3]-tstarts[3]
  
    pvalb = meb$cis$min.pv.gene
    pvalb = data.frame(names(pvalb), pvalb, stringsAsFactors=FALSE)  
    gns = sort(unique(pvalb[,1]))
    ords = pvalb
  
    nperm = ini.perm
    permind = 1:nperm
    permm = matrix(NA, nrow=nperm, ncol=nrow(pvalb))
    #do 1k perm
    starts = proc.time()
    for(i in permind){
      permi = sample(1:nsub)
      exprp = redboot[,permi]
      covap = permEst$cvrtM[permi,]
    
      cvrtp = SlicedData$new()
      cvrtp = cvrtp$CreateFromMatrix(t(covap))
        
      genep = SlicedData$new();
      genep = geneb$CreateFromMatrix(exprp)
      
      mep = Matrix_eQTL_main(
              snps = snps,
              gene = genep,
              cvrt = cvrtp,
              pvOutputThreshold = permEst$pvOutputThreshold,
              output_file_name = tmpbootfil,
              output_file_name.cis = bootfil,
              pvOutputThreshold.cis = permEst$pvOutputThreshold,
              useModel = useModel, 
              errorCovariance = errorCovariance,
              snpspos = permEst$snpspos,
              genepos = bootpos, 
              cisDist = permEst$cisDist,
              verbose = permEst$verbose,
              pvalue.hist = permEst$pvalue.hist,
              min.pv.by.genesnp = TRUE,
              noFDRsaveMemory = TRUE);
      pvalp = mep$cis$min.pv.gene
      pvalp = data.frame(names(pvalp), pvalp, stringsAsFactors=FALSE)  
      permm[i,] = pvalp[,2]
      message(i, " out of ", nperm)
    }
    ends = proc.time()
    tim0 = (ends[3]-starts[3])
    message("time per iter: ", tim0/nperm)#half minute per iteration per 1 gene
    colnames(permm) = pvalb[,1]
  
    y = rep(NA, ncol(permm))
    for(i in 1:ncol(permm)){
      y[i] = sum(permm[1:nperm,i]< pvalb[i,2])
    }
    message("attempt ", j, ": ", mean(y<(nperm*0.001)), " and ", mean(y>(nperm*.3)))
    tmpres[j, 1] = sum(y<(nperm*0.001))
    tmpres[j, 2] = sum(y>(nperm*.3))
    tmpres[j, 3] = tim0
  
  
    kp = (y/nperm>=0.001 & y/nperm<.3)
    kpboot = redboot[kp,]
    kpbootpos = bootpos[kp,]
    tmpres[j, 4] = sum(kp)
    if(sum(kp)>num.grid | j>=max.adj){
      break
    }
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
      }
    }
  message(j, " ", sum(kp))
  }
  for(recap in 1:j){
   message(tmpres[recap, 1], " ", tmpres[recap, 2], " ", tmpres[recap, 4], " ", tmpres[recap, 3], " ")
  }
  #
  #the rest of iterations
  #
  nperm0 = ini.perm+1
  nperm = n.perm
  if(nperm0<=nperm){
    permind = nperm0:nperm
    
    permm0 = permm
    permm = matrix(NA, nrow=nperm, ncol=nrow(pvalb))
    permm[1:ini.perm,] = permm0
  
  
    starts = proc.time()
    for(i in permind){
      permi = sample(1:nsub)
      exprp = redboot[,permi]
      covap = permEst$cvrtM[permi,]
  
      cvrtp = SlicedData$new()
      cvrtp = cvrtp$CreateFromMatrix(t(covap))
        
      genep = SlicedData$new();
      genep = geneb$CreateFromMatrix(exprp)
      
      mep = Matrix_eQTL_main(
              snps = snps,
              gene = genep,
              cvrt = cvrtp,
              pvOutputThreshold = permEst$pvOutputThreshold,
              output_file_name = tmpbootfil,
              output_file_name.cis = bootfil,
              pvOutputThreshold.cis = permEst$pvOutputThreshold,
              useModel = useModel, 
              errorCovariance = errorCovariance,
              snpspos = permEst$snpspos,
              genepos = bootpos, 
              cisDist = permEst$cisDist,
              verbose = permEst$verbose,
              pvalue.hist = permEst$pvalue.hist,
              min.pv.by.genesnp = TRUE,
              noFDRsaveMemory = TRUE);
      pvalp = mep$cis$min.pv.gene
      pvalp = data.frame(names(pvalp), pvalp, stringsAsFactors=FALSE)  
      permm[i,] = pvalp[,2]
      message(i, " out of ", nperm)
    }
    ends = proc.time()
    tim1 = (ends[3]-starts[3])
  }else{tim1=0}
  message("time per iter: ", (tim0+tim1)/nperm)#half minute per iteration per 1 gene



  y = rep(NA, ncol(permm))
  for(i in 1:ncol(permm)){
    y[i] = sum(permm[1:n.perm,i]< pvalb[i,2])
  }
  table(y)
  message("attempt ", j, ": ", mean(y<(nperm*0.001)), " and ", mean(y>nperm*.3))
  
  
  kp3 = (y/n.perm)>=0     & (y/n.perm)<=0.3
  kp3a = (y/n.perm)>0     & (y/n.perm)<=0.3
      
  y1 = -log10(y/nperm)
  x1 = -log10(pvalb[,2])
  glmi3 = glm(cbind(y[kp3],nperm-y[kp3])~x1[kp3], family="binomial")
  lmi3 = lm(y1[kp3a]~x1[kp3a])
  
  
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
  eigenMT$TESTSupd = round(pred.perm/eigenMT$p.value)
	cat(sprintf("pred.perm = %s\n",pred.perm))
	cat(sprintf("eigenMT$p.value = %s\n",eigenMT$p.value))
  if(eigenMT$TESTSupd<1)eigenMT$TESTSupd = 1
  if(eigenMT$TESTSupd>nrow(permEst$snpM))eigenMT$TESTSupd=nrow(permEst$snpM)
  
  if(!is.na(permEst$outdir)){
    updNtests = sprintf("%s/%s_updtests.csv", permEst$outdir, rownames(permEst$geneM)[1])
    filout1 = sprintf("%s/%s_boot_pval.csv", permEst$outdir, rownames(permEst$geneM)[1])
    filout2 = sprintf("%s/%s_short_boot_pval.csv", permEst$outdir, rownames(permEst$geneM)[1])
    filout3 = sprintf("%s/%s_time.csv", permEst$outdir, rownames(permEst$geneM)[1])
  }else{
    updNtests = sprintf("%s_updtests.csv", rownames(permEst$geneM)[1])
    filout1 = sprintf("%s_boot_pval.csv", rownames(permEst$geneM)[1])
    filout2 = sprintf("%s_short_boot_pval.csv", rownames(permEst$geneM)[1])
    filout3 = sprintf("%s_time.csv", rownames(permEst$geneM)[1])
  }
  if(write.res)
   write.csv(eigenMT, updNtests, quote=F, row.names=F)
  permm2 = rbind(permm, pvalb[,2])
  #ords$pvals,
  #ords$betas,
  #ords$tstat)
  colnames(permm2) = ords[,1]
  
  permp = rep(NA, ncol(permm))
  for(i in 1:ncol(permm)){
    permp[i] = mean(permm2[nperm+1,i]>=permm[1:nperm,i])
  }
  
  if(write.res){
    write.csv(permm2, filout1, row.names=F, quote=F) 
    write.csv(cbind(pvalb, permp), filout2, row.names=F, quote=F)
    write.csv(tim0+tim1, filout3, row.names=F, quote=F)
  }
  
  res = list(summ=eigenMT, tim=tim0+tim1, vals=cbind(pvalb, permp), mEQTL=me, min.snp = permEst$snpM[m,,drop=F])
  res
}
