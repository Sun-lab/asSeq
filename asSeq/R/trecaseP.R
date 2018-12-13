trecaseP <-
function(Y, Y1, Y2, X, Z, offset=NULL, 
      min.AS.reads=5, min.AS.sample=5, min.n.het=5, local.only=TRUE, 
      local.distance=200000, eChr=NULL, ePos=NULL, mChr=NULL, 
      mPos=NULL, converge=5e-5, convergeGLM=1e-8, scoreTestP=0.05, 
      transTestP=0.05, np.max=5000, np=c(20, 100, 500, 1000, 2500),
      aim.p=c(0.5, 0.2, 0.1, 0.05, 0.02), 
      confidence.p=0.01, trace=0, maxit=100)
{

  minVar = 1e-8
  
  if(!local.only){
    stop("sorry, local.only can only be TRUE\n")
  }
  
  if(min.AS.sample < min.n.het){
    stop("min.AS.sample should not be smaller than min.n.het\n")
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
  
  if(any(is.na(Y1))){
    stop("NA values in Y1\n")
  }
  
  if(any(is.na(Y2))){
    stop("NA values in Y2\n")
  }
  
  ## 
  # check Y1
  #
  
  if(is.vector(Y1)){
    nY = 1
    N1 = length(Y1)
  }else{
    nY = ncol(Y1)
    N1 = nrow(Y1)
  }
  
  ## 
  # check Y2
  #
  
  if(is.vector(Y2)){
    if(nY != 1) stop("dimensions of Y1 and Y2 do not match\n")
    
    N1 = length(Y2)
  }else{
    if(ncol(Y2) != nY) stop("dimensions of Y1 and Y2 do not match\n")
    if(nrow(Y2) != N1) stop("dimensions of Y1 and Y2 do not match\n")
  }
  
  ## 
  # check X, Y, and Z
  #
  
  if(is.vector(Y)){
    nY = 1
    N1 = length(Y)
    Y  = matrix(Y, ncol=1)
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
  
  if(! all(as.numeric(Z) %in% c(0,1,3,4))){
    stop("Z must take values 0, 1, 3, or 4\n")
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
  
  dims = numeric(9)
  dims[1] = nY
  dims[2] = nX
  dims[3] = nZ
  dims[4] = N1
  dims[5] = maxit
  dims[6] = useOffset
  dims[7] = min.AS.reads
  dims[8] = min.AS.sample
  dims[9] = min.n.het
  
  succeed = 0
    
  # First fit a baseline model using only the confouding covariates 
  # yFailBaselineModel indicate whether the baseline model can be fitted
  
  yFailBaselineModel = numeric(nY)
  
  # add one more column to X, whih is the space to be for one Z
  # in the glm computation
  
  X = cbind(X, rep(0,nrow(X)))
  
  # construct a genotype data matrix from the haplotype data matrix
  
  Zh = Z
  Z[Zh==3] = 1
  Z[Zh==4] = 2
  
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
  # best associated marker, and pvalue for TReC, ASE and TReCASE models 
  # permutation pvalue for TReC, ASE, TReCASE, and 
  # TReCASE or TReC depends on whether there it is local-eQTL or not
  
  bestM = rep(-9, 4*nY)
  pval  = rep(9,  4*nY) 
  perP  = numeric(4*nY)
  
  message("--------------------------------------------------------------")
  st1 = sprintf(" nY = %d, nX=%d, nZ=%d", ncol(Y), ncol(X)-1, ncol(Z))
  st1 = sprintf("%s, N=%d\n # of permutation=%d", st1, nrow(Y), nPermute)
  st1 = sprintf("%s, p-value threshhold=%.4e", st1, aim.p[1])
  message(st1)
  message("--------------------------------------------------------------")
  
  z1 = numeric(N1) 

  succeed = 0
  W1 = .C("trecase_permute", as.integer(dims), as.double(Y), as.double(X), 
          as.double(Z), as.double(z1), as.double(Y1), as.double(Y2), 
          as.double(Zh), as.integer(nPermute), as.double(offset), 
          as.integer(local.only), as.integer(local.distance), as.integer(eChr), 
          as.integer(ePos), as.integer(mChr), as.integer(mPos), 
          as.double(converge), as.double(convergeGLM), 
          yFailBaselineModel = as.integer(yFailBaselineModel), 
          as.double(scoreTestP), as.double(transTestP), 
          bestM = as.integer(bestM), pval=as.double(pval), 
          perP=as.double(perP), as.integer(trace), 
          succeed=as.integer(succeed), PACKAGE="asSeq")

  if(W1[["succeed"]] != 1){ stop("error in trecase_permute\n") }

  bestM = W1[["bestM"]]
  pval  = W1[["pval"]]
  perP  = W1[["perP"]]
  
  rm(W1)
  
  w3 = which(bestM < 0 | pval  < 0)
  if(length(w3) > 0) perP[w3] = NA
  
  npuse  = rep(nnp[1], length(perP))
  npuse[which(is.na(perP))] = NA

  bestM = matrix(bestM, ncol=4)
  pval  = matrix(pval,  ncol=4)
  perP  = matrix(perP,  ncol=4)  
  npuse = matrix(npuse, ncol=4)  
  
  # ---------------------------------------------------------
  # continue the permutations.
  # note that the last elements of np is the maximum number
  # of permutations to do. length(pcuts) = length(np) - 1
  # ---------------------------------------------------------
  
  for(i in 2:length(np)){
    which.kp = which(apply(perP, 1, min, na.rm=TRUE) < pcuts[i-1])
    
    if(length(which.kp)==0) break

    Ynew  =  Y[,which.kp,drop=FALSE]
    Y1new = Y1[,which.kp,drop=FALSE]
    Y2new = Y2[,which.kp,drop=FALSE]

    # update the number of genes
    dims[1] = ncol(Ynew)

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
    
    bbm  = rep(-9, 4*dims[1])
    pvs  = rep(9,  4*dims[1]) 
    ppvs = numeric(4*dims[1])
    
    message("--------------------------------------------------------------")
    st1 = sprintf(" nY=%d, nX=%d, nZ=%d", ncol(Ynew), ncol(X)-1, ncol(Z))
    st1 = sprintf("%s, N=%d\n # of permutation=%d", st1, nrow(Y), nPermute)
    if(i < length(np)){
      st1 = sprintf("%s, p-value threshhold=%.4e", st1, aim.p[i])
    }
    message(st1)
    message("--------------------------------------------------------------")

    succeed = 0
    z  = numeric(N1) 

    W2 = .C("trecase_permute", as.integer(dims), as.double(Ynew), as.double(X), 
           as.double(Z), as.double(z), as.double(Y1new), as.double(Y2new), 
           as.double(Zh), as.integer(nPermute), as.double(offset), 
           as.integer(local.only), as.integer(local.distance), as.integer(eChr1), 
           as.integer(ePos1), as.integer(mChr), as.integer(mPos), 
           as.double(converge), as.double(convergeGLM), 
           yFailBaselineModel = as.integer(yFailBaselineModel), 
           as.double(scoreTestP), as.double(transTestP), bestM = as.integer(bbm), 
           pval=as.double(pvs), perP=as.double(ppvs), as.integer(trace), 
           succeed=as.integer(succeed), PACKAGE="asSeq")
    
    if(W2[["succeed"]] != 1){ stop("error in trecase_permute\n") }
    
    bbm  = W2[["bestM"]]
    pvs  = W2[["pval"]]
    ppvs = W2[["perP"]]
    
    rm(W2)
    
    w3 = which(bbm < 0 | pvs < 0)
    if(length(w3) > 0) ppvs[w3] = NA
    
    bbm  = matrix(bbm,  ncol=4)
    pvs  = matrix(pvs,  ncol=4)
    ppvs = matrix(ppvs, ncol=4)  
    
    # ---------------------------------------------------------
    # update permutation p-values
    # w1 and w2 are weights for the previous permutations and  
    # the current permutations
    # ---------------------------------------------------------

    if(any(bbm != bestM[which.kp,,drop=FALSE])){
      stop("best associated markers do not match\n")
    }
    
    if(max(abs(as.numeric(pvs - pval[which.kp,,drop=FALSE])), na.rm=TRUE) > 1e-3){
      stop("best associated markers have different p-values\n")
    }
    
    for(j in 1:nrow(ppvs)){
      wj = which.kp[j]
      
      for(k in 1:ncol(ppvs)){
        
        w2 = nPermute/(nPermute + npuse[wj,k])
        w1 = 1 - w2

        if(is.na(ppvs[j,k]) && !(is.na(perP[wj,k]))){ 
          perP[wj,k]  = w1*perP[wj,k]
        }else if(is.na(perP[wj,k]) && !(is.na(ppvs[j,k]))){
          perP[wj,k]  = w2*ppvs[j,k]
        }else{
          perP[wj,k]  = w1*perP[wj,k] + w2*ppvs[j,k]
        }

        npuse[wj,k] = npuse[wj,k] + nPermute

      }
    }
    
  }

  gID   = 1:nY
  bestM[bestM < 0] = NA
  bestM = bestM + 1

  dataF = data.frame(geneID=gID, bestM, pval, perP, npuse)
  
  types = c("trec", "ase", "trecase", "trec_trecase")
  nms   = c("geneID", paste("markerID", types, sep="_"))
  nms   = c(nms, paste("pval", types, sep="_"))
  nms   = c(nms, paste("perP", types, sep="_"))
  nms   = c(nms, paste("nuse", types, sep="_"))
  
  names(dataF) = nms
  
  dataF

}

