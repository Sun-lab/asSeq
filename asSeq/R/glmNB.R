glmNB <-
function(y, X, link="log", offset=NULL, naPercent=0.4, maxit=50,
          conv=1e-8, fitted=NULL, scoreTestP=0.05, trace=1)
{
  if(!is.numeric(y)){
    stop("y must be a numeric vector\n")
  }
  
  if(!is.matrix(X)){
    stop("X must be a vector\n")
  }
  
  M = ncol(X)
  N = length(y)
  if(nrow(X) != N){
    stop("the dimension of X and y do not match\n")
  }
  
  if(!is.null(offset)){
    useOffset = 1
    if((! is.numeric(offset)) || length(offset) != N){
      stop("offset must be a numeric vector of the same length as y\n")
    }
  }else{
    useOffset = 0
    offset    = 0.0 
  }
    
  if(is.null(fitted)){
    init = 0
    fitted = rep(0, N)
  }else{
    init = 1

    if((! is.numeric(fitted)) || length(fitted) != N){
      stop("fitted must be a numeric vector of the same length as y\n")
    }
  }
    
  isNA = apply(X, 1, function(v){ any(is.na(v)) })
  isNA = is.na(y) | isNA
  
  if(length(which(isNA))/N > naPercent){
    stop("percent of missing data is too high\n")
  }
  
  w2kp  = which(!isNA)
  yR    = y[w2kp]
  XR    = X[w2kp,]
  
  if(useOffset){
    offset = offset[w2kp]
  }

  N  = length(yR)
  
  dims = numeric(5)
  dims[1] = N
  dims[2] = M
  dims[3] = maxit
  dims[4] = init
  dims[5] = useOffset
  
  rank    = scale = 1
  Xb      = XR
  resid   = weights = numeric(N)
  fitted  = fitted[w2kp]
  dfResid = 0
  nIter   = 0
  phi     = 0.0
  beta    = 0.0
  twologlik  = 0.0

  # NB Family
  fam0 = 5

  if(link=="logit"){
    linkR = 1
  }else if(link=="log"){
    linkR = 2
  }else if(link=="identity"){
    linkR = 3
  }else if(link=="inverse"){
    linkR = 4
  }else{
    stop("invalid family\n")
  }
  
  # z is working version of y
  z = numeric(N) 

  Z = .C("glmNB", as.integer(dims), nIter=as.integer(nIter), as.double(yR), 
         as.double(z), as.integer(linkR), offset=as.double(offset), 
         as.double(XR), as.double(conv), rank=as.integer(rank), as.double(Xb), 
         fitted=as.double(fitted), resid=as.double(resid), 
         weights=as.double(weights), phi=as.double(phi),
         scale=as.double(scale), dfResid=as.integer(dfResid), 
         as.integer(fam0), twologlik=as.double(twologlik), 
         as.double(scoreTestP), as.integer(trace), as.double(beta), 
         PACKAGE="asSeq")
  
  N = length(y)
  fitted = weights = resid = rep(NA, N)
  fitted[w2kp]  = Z$fitted
  weights[w2kp] = Z$weights
  resid[w2kp]   = Z$resid
  
  if(useOffset){
    offset = rep(NA, N)
    offset[w2kp] = Z$offset
  }

  xx = list(phi=Z$phi, offset=offset, loglik=Z$twologlik/2,
            rank=Z$rank, fitted=fitted, resid=resid, weights=weights, 
            scale=Z$scale, dfResid=Z$dfResid, nIter=Z$nIter, X=X)
  
  class(xx) = "glmNB"
  xx
  
}

