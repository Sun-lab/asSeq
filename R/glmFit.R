glmFit <-
function(family, link, y,  X, offset=NULL, naPercent=0.4, nTotal=NULL, 
         maxit=50, conv=1e-5, fitted=NULL, trace=1)
{
  if(!is.numeric(y)){
    stop("y must be a numeric vector\n")
  }
  
  if(!is.matrix(X)){
    stop("X must be a matrix\n")
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
  
  if(family == "binomial"){
    if((! is.numeric(nTotal)) || length(nTotal) != N){
      strWarn = "For binomial family, nTotal must be a numeric vector"
      stop(strWarn, "of the same length as y\n")
    }
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
  
  #define BINOMIAL  1
  #define POISSON   2
  #define GAUSSIAN  3
  #define GAMMA     4
  
  if(family=="binomial"){
    familyR = 1
  }else if(family=="poisson"){
    familyR = 2
  }else if(family=="gaussian"){
    familyR = 3
  }else if(family=="gamma"){
    familyR = 4
  }else{
    stop("invalid family\n")
  }
  
  #define LOGIT     1
  #define LOG       2
  #define IDENTITY  3
  #define INVERSE   4
  
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
  
  isNA = apply(X, 1, function(v){ any(is.na(v)) })
  isNA = is.na(y) | isNA
  
  if(length(which(isNA))/N > naPercent){
    stop("percent of missing data is too high\n")
  }
  
  w2kp = which(!isNA)
  yR = y[w2kp]
  XR = X[w2kp,]
  
  if(useOffset){
    offset = offset[w2kp]
  }
  if(family == "binomial"){
    nTotal = nTotal[w2kp]
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
  phi     = 0
  beta    = 995.0
  # z is working version of y
  z       = numeric(N) 
  
  Z = .C("glmFit", as.integer(familyR), as.integer(linkR), 
         as.integer(dims), nIter=as.integer(nIter), 
         as.double(yR), offset=as.double(offset), as.double(z),
         as.double(XR), as.double(nTotal), as.double(conv), 
         rank=as.integer(rank), as.double(Xb), fitted=as.double(fitted), 
         resid=as.double(resid), weights=as.double(weights), 
         as.double(phi), as.integer(trace), scale=as.double(scale), 
         dfResid=as.integer(dfResid), beta=as.double(beta), PACKAGE="asSeq")
  
  N = length(y)
  fitted = weights = resid = rep(NA, N)
  fitted[w2kp]  = Z$fitted
  weights[w2kp] = Z$weights
  resid[w2kp]   = Z$resid
  
  if(useOffset){
    offset = rep(NA, N)
    offset[w2kp] = Z$offset
  }

  xx = list(faimily=family, link=link, offset=offset, 
            rank=Z$rank, fitted=fitted, resid=resid, weights=weights, 
            scale=Z$scale, dfResid=Z$dfResid, nIter=Z$nIter, X=X, 
            beta=Z$beta)
  
  class(xx) = "glmFit"
  xx
  
}

