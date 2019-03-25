ase <-
  function(Y1, Y2, Z, geneloc, SNPloc, file_ase="ase.txt", cis_window = 1e5L,
           min_ASE_total = 8L, min_nASE = 10L,
           eps = 5e-5, max_iter = 400L, show = FALSE)
  {
    
    minVar = 1e-8
    
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
    
    if(! all(as.numeric(Z) %in% c(0,1,2,3))){
      stop("Z must take values 0, 1, 2, or 3\n")
    }
    
    ## ----------------------------
    ## check the SNP and gene location 
    ## ----------------------------
    
    if(nrow(geneloc) != nY)
      stop("number of Gene in GeneInfo doesn't match that in Y")
    if(nrow(SNPloc) != nZ)
      stop("number of SNP in GeneInfo doesn't match that in Z")
    
    if(!all(colnames(SNPloc) == c("snp", "chr", "pos"))){
      stop("colnames of SNPloc has to be c('snp', chr, 'pos')")
    }
    
    if(!all(colnames(geneloc) == c("gene", "chr", "start", 'end'))){
      stop("colnames of SNPloc has to be c('gene', chr, 'start','end')")
    }
  
    # sort geneloc and snploc by chr and pos 
    geneloc$chr = gsub("chr", "", geneloc$chr)
    geneloc$chr[which(geneloc$chr == "X")] = 23
    geneloc$chr[which(geneloc$chr == "Y")] = 24
    
    geneloc$chr   = as.integer(geneloc$chr)
    geneloc$start = as.integer(geneloc$star)
    geneloc$end   = as.integer(geneloc$end)
    
    SNPloc$chr = gsub("chr", "", SNPloc$chr)
    SNPloc$chr[which(SNPloc$chr == "X")] = 23
    SNPloc$chr[which(SNPloc$chr == "Y")] = 24
    SNPloc$chr = as.integer(SNPloc$chr)
    SNPloc$pos = as.integer(SNPloc$pos)
    
    
    # check chr number are at same order 
    if(!all(unique(SNPloc$chr) == unique(geneloc$chr))){
      stop("SNPloc and geneloc does not have same chromosome order.")
    }
    
    Rcpp_ase_mtest(Y1, Y2, ZZ, SNPloc$pos, SNPloc$chr, 
                   geneloc$start, geneloc$end, geneloc$chr, 
                   file_ase, cis_window, min_ASE_total, 
                   min_nASE, eps, max_iter, show)
    
  }
