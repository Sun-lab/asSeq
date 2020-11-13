get_reduced_boot = function(ppi,target.min.ps=NULL,target.perm.ps, genenm=NULL, mQTL.fit, expr.mat, min.SNP, covars, nsam, use.norm=F){
  if(is.null(genenm))genenm=rownames(expr.mat)[1]
  if(is.null(genenm))genenm=1:length(genenm)
  #ppi = 50;i = 1; mQTL.fit=eigenMT; expr.mat = exprj; min.SNP=gen.sub; target.perm.p=target.perm.pi;use.norm=F; covars = covar
  #i = i + 1
  i = 1
  target.perm.p = target.perm.ps[ppi]
  #nsam = ncol(gen.sub)-1
  if(is.null(target.min.ps)){
    target.min.p = target.perm.p/mQTL.fit$TESTS[i]#/mQTL.fit$ntest[i]#/mQTL.fit$TESTS[i]
  }else{
    target.min.p = target.min.ps[ppi]
  }
  tstat = qt(target.min.p/2, df=nsam-ncol(covars)-1, lower.tail=F)
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
	
	# These three lines aren't used elsewhere 
	#		in this function and are throwing an 
	#		error when the min.SNP contains NAs ...
  # ypre0 = fitspre+lm1$residuals
  # d2 = data.frame(y=ypre0, covars, xeqtl)
  # lm2pre0 = lm(y~., data=d2)
	
	# If there are missing genotypes, these 
	#		lines will throw an error
	idx_nonNA = !is.na(xeqtl)
	nn = nrow(covars)
	ypre = rep(NA,nn)
	ypre[idx_nonNA] = fitspre+res
  # ypre = fitspre+res
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
  rownames(yout) = sprintf("%s_%s", genenm, ppi)
  #summary(lm2out)
  #target.min.p
#  colnames(yout) = sprintf("%s_%s", rownames(expr.mat)[i], ppi)
  yout
}

