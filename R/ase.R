ase <-
function(Y1, Y2, Z, output.tag, p.cut, min.AS.reads=5, min.AS.sample=5, 
  min.n.het=5, local.only=TRUE,  local.distance=2e+05, 
  eChr=NULL, ePos=NULL, mChr=NULL, mPos=NULL, trace=0)
{
  
  minVar = 1e-8
  
  if(min.AS.reads < min.n.het){
    stop("min.AS.reads should not be smaller than min.n.het\n")
  }
  
  ## 
  # check the NAs in Y and Z
  #
  if(any(is.na(Y1))){
    stop("NA values in Y1\n")
  }
  
  if(any(is.na(Y2))){
    stop("NA values in Y2\n")
  }
  
  if(any(is.na(Z))){
    stop("NA values in Z\n")
  }
  
  ## 
  # check Y1
  #
  
  if(is.vector(Y1)){
    nY = 1
    N1 = length(Y1)
    if(var(Y1) < minVar){
      stop("Variance of Y1 is close to 0\n")
    }
  }else{
    nY = ncol(Y1)
    N1 = nrow(Y1)
    
    varY = apply(Y1, 2, var)
    wVar = which(varY < minVar)
    
    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Y1 have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
  }
  
  ## 
  # check Y2
  #
    
  if(is.vector(Y2)){
    if(nY != 1) stop("dimensions of Y1 and Y2 do not match\n")
    
    N1 = length(Y2)
    if(var(Y2) < minVar){
      stop("Variance of Y2 is close to 0\n")
    }
  }else{
    if(ncol(Y2) != nY) stop("dimensions of Y1 and Y2 do not match\n")

    if(nrow(Y2) != N1) stop("dimensions of Y1 and Y2 do not match\n")
    
    varY = apply(Y2, 2, var)
    wVar = which(varY < minVar)
    
    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Y2 have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
  }

  ## 
  # check Z
  #
  
  if(is.vector(Z)){
    nZ = 1
    N2 = length(Z)
    
    if(var(Z) < 1e-7){
      stop("Variance of X is 0\n")
    }
  }else{
    nZ = ncol(Z)
    N2 = nrow(Z)
    
    varZ = apply(Z, 2, var)
    wVar = which(varZ < minVar)
    
    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Z have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
  }
  
  if(N1 != N2){
    stop("X and Z must have the same number of rows!")
  }
  
  if(! all(as.numeric(Z) %in% c(0,1,3,4))){
    stop("Z must take values 0, 1, 3, or 4\n")
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
  }else{
    warning("ASE mapping requires haplotype information, which may be inaccruate for non-local eQTL\n")
  }
  
  dims = numeric(6)
  dims[1] = nY
  dims[2] = nZ
  dims[3] = N1
  dims[4] = min.AS.sample
  dims[5] = min.AS.reads
  dims[6] = min.n.het

  output = character(2)
  output[1] = sprintf("%s_eqtl.txt", output.tag)
  output[2] = sprintf("%s_freq.txt", output.tag)
    
  # Now, go for it
  succeed = 0
  
  W = .C("ase", as.integer(dims), as.double(Y1), as.double(Y2), 
         as.double(Z), as.character(output), as.double(p.cut), 
         as.integer(local.only), as.integer(local.distance), 
         as.integer(eChr), as.integer(ePos), 
         as.integer(mChr), as.integer(mPos), 
         as.integer(trace), succeed=as.integer(succeed), 
         PACKAGE="asSeq")
  
  list(succeed = W[["succeed"]])
  
}

