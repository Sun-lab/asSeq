trec <-
function(Y, X, Z, output.tag, p.cut, link="log", offset=NULL, 
      adjZ=TRUE, local.only=TRUE, local.distance=2e+05, eChr=NULL, ePos=NULL,
      mChr=NULL, mPos=NULL, converge=5e-5, convergeGLM=1e-8, scoreTestP=0.05, 
      trace=1, maxit=100)
{
  
  if(link=="log"){
    linkR = 2  
  }else if(link=="identity"){
    linkR = 3
    if(adjZ){
      adjZ = FALSE
      warning("adjZ is set to be FALSE for identity link\n")
    }
  }else{
    stop("only 'log' or 'identity' link functions are allowed\n")
  }
  
  ## 
  # check the NAs in X, Y, and Z
  #
  if(any(is.na(Y))){
    stop("NA values in Y\n")
  }
  
  if(any(is.na(X))){
    stop("NA values in X\n")
  }
  
  if(any(is.na(Z))){
    stop("NA values in Z\n")
  }
  
  ## 
  # check the dimensions of X, Y, and Z
  #
  if(is.vector(Y)){
    nY = 1
    N1 = length(Y)
    if(var(Y) < converge){
      stop("Variance of Y is close to 0\n")
    }
  }else{
    nY = ncol(Y)
    N1 = nrow(Y)
    
    varY = apply(Y, 2, var)
    wVar = which(varY < converge)
    
    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Y have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
  }
  
  if(is.vector(X)){
    nX = 1
    N2 = length(X)

    if(var(X) < converge){
      stop("Variance of X is 0\n")
    }
  }else{
    nX = ncol(X)
    N2 = nrow(X)
    
    varX = apply(X, 2, var)
    wVar = which(varX < converge)
    
    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in X have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
  }
  
  if(is.vector(Z)){
    nZ = 1
    N3 = length(Z)
    
    if(var(Z) < converge){
      stop("Variance of X is 0\n")
    }
  }else{
    nZ = ncol(Z)
    N3 = nrow(Z)
    
    varZ = apply(Z, 2, var)
    wVar = which(varZ < converge)
    
    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Z have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
  }
    
  if(N1 != N2 || N1 != N3){
    stop("X, Y, and Z must have the same number of rows!")
  }
  
  ## 
  # check the length of offset
  #
  if(!is.null(offset)){
    if(length(offset) != N1){
      stop("length of offset does not match the sample size in Y\n")
    }
    useOffset = 1
  }else{
    useOffset = 0
    offset = rep(0, N1)
  }
  
  ## 
  # check the length of chromosome location information
  #
  if(local.only){
    if(!is.numeric(eChr) || length(eChr) != nY){
      stop(sprintf("eChr must be a numeric vector of length %d\n", nY))
    }
    if(!is.numeric(ePos) || length(ePos) != nY){
      stop(sprintf("ePos must be a numeric vector of length %d\n", nY))
    }
    if(!is.numeric(mChr) || length(mChr) != nZ){
      stop(sprintf("mChr must be a numeric vector of length %d\n", nZ))
    }
    if(!is.numeric(mPos) || length(mPos) != nZ){
      stop(sprintf("mPos must be a numeric vector of length %d\n", nZ))
    }
  }
  
  dims = numeric(6)
  dims[1] = nY
  dims[2] = nX
  dims[3] = nZ
  dims[4] = N1
  dims[5] = maxit
  dims[6] = useOffset
  
  succeed = 0
  
  output = character(2)
  output[1] = sprintf("%s_eqtl.txt", output.tag)
  output[2] = sprintf("%s_freq.txt", output.tag)

  # First fit a baseline model using only the confouding covariates 
  # yFailBaselineModel indicate whether the baseline model can be fitted
  
  yFailBaselineModel = numeric(nY)
  
  # add one more column to X, whih is the space to be for one Z
  # in the glm computation
  
  X = cbind(X, rep(0,nrow(X)))
  # z is working version of y
  z1 = numeric(N1) 

  W = .C("glmEQTL", as.integer(dims), as.double(Y), as.double(X), as.double(Z), 
         as.double(z1), as.integer(linkR), as.double(offset), 
         as.integer(adjZ), as.character(output), as.double(p.cut), 
         as.integer(local.only), as.integer(local.distance), 
         as.integer(eChr), as.integer(ePos), as.integer(mChr), 
         as.integer(mPos), as.double(converge), as.double(convergeGLM), 
         yFailBaselineModel = as.integer(yFailBaselineModel), 
         as.double(scoreTestP), as.integer(trace), 
         succeed=as.integer(succeed), PACKAGE="asSeq")
    
  list(succeed = W[["succeed"]], yFailBaselineModel=W[["yFailBaselineModel"]])
  
}

