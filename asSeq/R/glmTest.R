glmTest <-
function(gm, Z)

{
  if(class(gm) != "glmFit"){
    stop("gm must be of class glmFit\n")
  }
  
  resid   = gm$resid
  weights = gm$weights
  scale   = gm$scale
  
  N = length(resid)
  
  isNA1 = is.na(resid)
  isNA2 = is.na(weights)
  
  if(any(isNA1 != isNA2)){
    stop("NAs in ressid and weights do not match")
  }
  
  if(!is.matrix(Z)){
    stop("Z must be a vector\n")
  }
  
  P = ncol(Z)
  if(nrow(Z) != N){
    stop("the dimension of Z does not match the dimension of the model\n")
  }
  
  w2kp  = which(!isNA1)
  XR    = gm$X[w2kp,,drop=FALSE]
  resid = resid[w2kp]
  weights = weights[w2kp]
  
  isNA = apply(Z, 1, function(v){ any(is.na(v)) })
  if(any(isNA)){ stop("missing value in Z\n") }
  
  N = length(resid)
  M = ncol(XR)
  P = ncol(Z)
  
  dims = c(N, M, P)
  chi2 = 0.0
  df   = 0
  
  W = .C("glm_score_test", as.integer(dims), as.double(Z), 
         as.double(resid), as.double(weights), as.double(XR), 
         as.double(scale), chi2=as.double(chi2), 
         df = as.integer(df), PACKAGE="asSeq")
  
  chi2 = W$chi2
  df   = W$df
  pval = pchisq(chi2, df, lower.tail=FALSE)
  
  list(chi2=chi2, df=df, pval=pval)
}

