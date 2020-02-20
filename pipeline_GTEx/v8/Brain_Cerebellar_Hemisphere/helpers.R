get_blocks = function(x, split="\\|", blocks=c(1)){
  paste(unlist(strsplit(x, split=split))[blocks], collapse="-")
}

minord = function(genei, pvals, genes){
  kp = which(genes==genei)
  c(as.character(genei), kp[order(pvals[kp])[1]])
}

classify = function(i, frac1, fracN, kp){
  out = c(mean(abs(frac1[i,kp[i,]]-.5)>.4), mean(fracN[i,kp[i,]]), sum(kp[i,]))
}


get_reduced_boot = function(ppi,target.min.ps=NULL,target.perm.ps, i=1, mQTL.fit, expr.mat, min.SNP, covars, nsam, use.norm=F){
  #ppi = 50;i = 1; mQTL.fit=eigenMT; expr.mat = exprj; min.SNP=gen.sub; target.perm.p=target.perm.pi;use.norm=F; covars = covar
  #i = i + 1
  target.perm.p = target.perm.ps[ppi]
  nsam = ncol(gen.sub)-1
  if(is.null(target.min.ps)){
    target.min.p = target.perm.p/mQTL.fit$TESTS[i]#/mQTL.fit$ntest[i]#/mQTL.fit$TESTS[i]
  }else{
    target.min.p = target.min.ps[ppi]
  }
  tstat = qt(target.min.p/2, df=nsam-7, lower.tail=F)
  xeqtl = as.numeric(min.SNP[i,])
  #x1 = as.numeric(covars[1,])
  #x2 = as.numeric(covars[2,])
  #x3 = as.numeric(covars[3,])
  #x4 = as.numeric(covars[4,])
  #x5 = as.numeric(covars[5,])
  y = as.numeric(expr.mat[i,])
  d = data.frame(y=y, covars, xeqtl)
  lm1 = lm(y~., data=d)
  npar = length(lm1$coef)
  se = summary(lm1)$coefficients[npar,2]
  if(use.norm){
    sdi = sd(lm1$residuals)
    res = rnorm(nsam, 0, sdi)
  }else{
    res = sample(lm1$residuals)
  }
  fitspre = lm1$fitted
  ypre0 = fitspre+lm1$residuals
  d2 = data.frame(y=ypre0, covars, xeqtl)
  lm2pre0 = lm(y~., data=d2)

  ypre = fitspre+res
  d3 = data.frame(y=ypre, covars, xeqtl)
  lm2pre = lm(y~., data=d3)
  respre = lm2pre$residuals

  seout = summary(lm2pre)$coefficients[npar,2]
  cfout = lm2pre$coef
  cfout[npar] = seout*tstat#cf1["xeqtl"]*mult
  fitsout = c(cbind(rep(1,nrow(covars)),covars,xeqtl)%*%cfout)
  yout = matrix(fitsout + respre, nrow=1)

  d4 = data.frame(y=yout[1,], covars, xeqtl)
  lm2out = lm(y~., data=d4)
  rownames(yout) = sprintf("%s_%s", rownames(expr.mat)[i], ppi)
  #summary(lm2out)
  #target.min.p
#  colnames(yout) = sprintf("%s_%s", rownames(expr.mat)[i], ppi)
  yout
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

