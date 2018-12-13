trecP <-
function(Y, X, Z, link="log", weight=NULL, offset=NULL, adjZ=TRUE, 
  local.only=TRUE, local.distance=2e+05, eChr=NULL, ePos=NULL,
  mChr=NULL, mPos=NULL, converge=1e-8, scoreTestP=0.05, 
  np.max=5000, np=c(20, 100, 500, 1000, 2500),
  aim.p=c(0.5, 0.2, 0.1, 0.05, 0.02), 
  confidence.p=0.01, trace=0, maxit=100)
{
  if(link=="log"){
    linkR = 2  
  }else if(link=="identity"){
    linkR = 3
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
    
    if(var(X) < 1e-7){
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
    
    if(var(Z) < 1e-7){
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
  # check the length of prior weight and offset
  #
  if(!is.null(weight)){
    if(length(weight) != N1){
      stop("length of weight does not match the sample size in Y\n")
    }
  }else{
    weight = rep(1, N1)
  }
  
  if(!is.null(offset)){
    if(length(offset) != N1){
      stop("length of offset does not match the sample size in Y\n")
    }
    useOffset = 1
  }else{
    useOffset = 0
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
    
  # First fit a baseline model using only the confouding covariates 
  # yFailBaselineModel indicate whether the baseline model can be fitted
  
  yFailBaselineModel = numeric(nY)
  
  # add one more column to X, whih is the space to be for one Z
  # in the glm computation
  
  X = cbind(X, rep(0,nrow(X)))
  
  # ---------------------------------------------------------
  # for each gene, evaluate the pemrutation p-value of its
  # most significant nominal p-value
  # ---------------------------------------------------------

  if(length(np) != length(aim.p)){
    stop("the length of np must equal to the length of aim.p\n")
  }

  if(any(diff(np)<=0)){ stop("np must be strictly ascending\n") }

  if(any(diff(aim.p)>=0)){ stop("aim.p must be strictly descending\n") }

  if(np.max <= max(np)){
    stop("np.max must be bigger than max(np)\n")
  }

  cuts  = qbinom(1-confidence.p, np, aim.p)
  pcuts = cuts/np
  
  # np:  total number of permutations for each targeted p-value
  # nnp: number of additional permutations for each targeted p-value
  np  = c(np, np.max)
  nnp = c(np[1], diff(np))
        
  # ---------------------------------------------------------
  # first carry out "nnp[1]" permutations. Record the best 
  # associated markers and corresponding nominal/permuted 
  # p-value for eah gene.
  #
  # pval: nominal p-value
  # perP: permutation p-value
  # ---------------------------------------------------------
  
  # the number of permutatios to do
  nPermute = nnp[1]
  
  ###
  # best associated marker, pvalue and permutation pvalue of 
  # the best association
  
  bestM = rep(-9, nY)
  pval  = rep(-1,  nY) 
  perP  = numeric(nY)
  
  message("--------------------------------------------------------------")
  st1 = sprintf(" nY = %d, nX=%d, nZ=%d", ncol(Y), ncol(X), ncol(Z))
  st1 = sprintf("%s, N=%d\n # of permutation=%d", st1, nrow(Y), nPermute)
  st1 = sprintf("%s, p-value threshhold=%.4e", st1, aim.p[1])
  message(st1)
  message("--------------------------------------------------------------")

  W = .C("glmEQTL_permute", as.integer(dims), as.double(Y), as.double(X), 
         as.double(Z), as.integer(nPermute), as.integer(linkR), as.double(weight), 
         as.double(offset), as.integer(adjZ),  as.integer(local.only), 
         as.integer(local.distance), as.integer(eChr), as.integer(ePos), 
         as.integer(mChr), as.integer(mPos), as.double(converge), 
         yFailBaselineModel = as.integer(yFailBaselineModel), 
         as.double(scoreTestP), bestM = as.integer(bestM), 
         pval=as.double(pval), perP=as.double(perP), as.integer(trace), 
         succeed=as.integer(succeed), PACKAGE="asSeq")

  bestM = W[["bestM"]]
  pval  = W[["pval"]]
  perP  = W[["perP"]]
  
  w1 = which(bestM < 0)
  if(length(w1) > 0) bestM[w1] = NA
  
  w2 = which(pval < 0)
  if(length(w2) > 0) pval[w2] = NA

  w3 = union(w1, w2)
  if(length(w3) > 0) perP[w3] = NA
  
  # number of permutations for each gene
  npall  = rep(nnp[1], nY)
  
  # ---------------------------------------------------------
  # continue the permutations.
  # note that the last elements of np is the maximum number
  # of permutations to do. length(pcuts) = length(np) - 1
  # ---------------------------------------------------------
  
  for(i in 2:length(np)){
    which.kp = which(perP < pcuts[i-1])
    
    if(length(which.kp)==0) break

    Y1 = Y[,which.kp,drop=FALSE]
    
    # update the number of genes
    dims[1] = ncol(Y1)

    # update the number of permutations to do
    nPermute = nnp[i] 

    # update the chromosome information
    if(local.only){
      eChr1 = eChr[which.kp]
      ePos1 = ePos[which.kp]
    }else{
      eChr1 = 0
      ePos1 = 0
    }
    
    bbm = pvs = ppvs = numeric(dims[1])
    
    message("--------------------------------------------------------------")
    st1 = sprintf(" nY = %d, nX=%d, nZ=%d", ncol(Y1), ncol(X), ncol(Z))
    st1 = sprintf("%s, N=%d\n # of permutation=%d", st1, nrow(Y1), nPermute)
    if(i < length(np)){
      st1 = sprintf("%s, p-value threshhold=%.4e", st1, aim.p[i])
    }
    message(st1)
    message("--------------------------------------------------------------")

    W = .C("glmEQTL_permute", as.integer(dims), as.double(Y1), as.double(X), 
           as.double(Z), as.integer(nPermute), as.integer(linkR), as.double(weight), 
           as.double(offset), as.integer(adjZ), as.integer(local.only), 
           as.integer(local.distance), as.integer(eChr1), as.integer(ePos1), 
           as.integer(mChr), as.integer(mPos), as.double(converge), 
           yFailBaselineModel = as.integer(yFailBaselineModel), 
           as.double(scoreTestP), bestM = as.integer(bbm), 
           pval=as.double(pvs), perP=as.double(ppvs), as.integer(trace), 
           succeed=as.integer(succeed), PACKAGE="asSeq")
    
    # ---------------------------------------------------------
    # update permutation p-values
    # w1 and w2 are weights for the previous permutations and  
    # the current permutations
    # ---------------------------------------------------------

    if(any(W$bestM != bestM[which.kp])){
      stop("best associated markers do not match\n")
    }
    
    w2 = nnp[i]/np[i]
    w1 = 1 - w2
    perP[which.kp] = w1*perP[which.kp] + w2*W[["perP"]]
    npall[which.kp] = np[i]
  }

  gID = 1:nY
  mID = bestM + 1
  
  dataF = data.frame(geneID=gID, markerID=mID, pval=pval,
    permuteP=perP, npermute=npall)
  
  return(dataF)

}

