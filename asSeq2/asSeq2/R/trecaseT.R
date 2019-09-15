trecaseT <-
  function(Y, Y1 = NULL, Y2 = NULL, Z, XX, RHO, CNV1, CNV2,
           SNPloc, geneloc, GeneSnpList = list(),
           file_trec = "trecT.txt", file_trecase = "trecaseT.txt", 
           useLRT = FALSE, transTestP = 0.01, cis_window = 100000, useASE = 1L, 
           min_ASE_total = 8L, min_nASE = 5L, min_nASE_het = 5L, eps = 0.00001, 
           max_iter = 4000L, show = FALSE)
  {
    ## Y: matrix of total read count. Each row is a sample, and each column is a
    ##    gene
    ## Y1, Y2: matrix of allele-specific read count. Each row is a sample, and 
    ##    each column is a gene 
    ## Z: matrix of genotype data. Each row is a sample, and each column is a SNP
    ## XX: covariates needs to be adjusted in Negtive Bionomial Regression
    ## SNPloc/geneloc: data.frame of SNP location information, the column names 
    ##    have to be c("snp", "chr", "pos") / c("gene", "chr", "start","end")
    ## file_trec/file_trecase: output file name of trec/trecase 
    ## cis_window:
    ## useASE:
    ## min_ASE_total:
    ## min_nASE:
    ## eps:
    ## max_iter:
    ## show:
    
    minVar = 1e-8
    
    ## ----------------------------
    ## check the NAs 
    ## ----------------------------
    
    if(any(is.na(Y))){
      stop("NA values in Y\n")
    }
    
    if(any(is.na(XX))){
      stop("NA values in X\n")
    }
    if(any(is.na(RHO))){
      stop("NA values in RHO\n")
    }
    if(any(is.na(Z))){
      Z[is.na(Z)] = -9
      warning("NA values in Z\n")
    }
    if(any(is.na(CNV1))){
      CNV1[is.na(CNV1)] = -9
      warning("NA values in CNV1\n")
    }
    if(any(is.na(CNV2))){
      CNV2[is.na(CNV2)] = -9
      warning("NA values in CNV2\n")
    }
    
    ## ----------------------------
    ## change the format X, Y, and Z
    ## ----------------------------
    
    Y    = data.matrix(Y)
    Z    = data.matrix(Z)
    XX   = data.matrix(XX)
    CNV1 = data.matrix(CNV1)
    CNV2 = data.matrix(CNV2)
    
    nGene = ncol(Y)
    nSam  = nrow(Y)
    
    ## ----------------------------
    ## make sure the dims consistant
    ## ----------------------------
    
    if(nrow(Z) != nSam)
      stop("Z and Y have different sample size")
    
    if(nrow(CNV1) != nSam | nrow(CNV2) != nSam)
      stop("CNV and Y have different sample size")
    
    if(ncol(CNV1) != nGene | ncol(CNV2) != nGene)
      stop("CNV and Y have different number of genes")
    
    if(length(RHO) != nSam)
      stop("the length of RHO doesn't match sample size")
    
    ## ----------------------------
    ## check the variance 
    ## ----------------------------
    varY = apply(Y, 2, var)
    wVar = which(varY < minVar)
    
    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Y have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
    
    # add intercept if there isn't
    if(any(abs(XX[,1]-1) > 1e-8)){
      XX = model.matrix( ~ XX)  
    }
    
    # check the variance of XX    
    varX  = apply(XX, 2, var)
    wVarX = which(varX < minVar)
    if(any(wVarX > 1)){
      stop("XX has tiny variance")
    }
    
    varZ = apply(Z, 2, var)
    wVar = which(varZ < minVar)
    
    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Z have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
    
    if(! all(as.numeric(Z) %in% c(0,1,2,3, -9))){
      stop("Z must take values 0, 1, 2, 3 or NA\n")
    }
    
    ## ----------------------------
    ## check the SNP and gene location 
    ## ----------------------------
    
    if(nrow(geneloc) != nGene)
      stop("number of Gene in GeneInfo doesn't match that in Y")
    if(nrow(SNPloc) != ncol(Z))
      stop("number of SNP in GeneInfo doesn't match that in Z")
    
    if(!all(colnames(SNPloc) == c("snp", "chr", "pos"))){
      stop("colnames of SNPloc has to be c('snp', chr, 'pos')")
    }
    
    if(!all(colnames(geneloc) == c("gene", "chr", "start", 'end'))){
      stop("colnames of geneloc has to be c('gene', chr, 'start','end')")
    }
    
    # geneloc has to be sorted by chromosome
    geneloc$chr = gsub("chr", "", geneloc$chr)
    geneloc$chr[which(geneloc$chr == "X")] = 23
    geneloc$chr[which(geneloc$chr == "Y")] = 24
    
    geneloc$chr   = as.integer(as.character(geneloc$chr))
    geneloc$start = as.integer(as.character(geneloc$star))
    geneloc$end   = as.integer(as.character(geneloc$end))
    
    SNPloc$chr = gsub("chr", "", SNPloc$chr)
    SNPloc$chr[which(SNPloc$chr == "X")] = 23
    SNPloc$chr[which(SNPloc$chr == "Y")] = 24
    SNPloc$chr = as.integer(as.character(SNPloc$chr))
    SNPloc$pos = as.integer(as.character(SNPloc$pos))
    
    if(any(sort(geneloc$chr) != geneloc$chr)){
      stop("geneloc has to be sorted by chromosome\n") 
    }
    
    if(useASE){
      
      ## ----------------------------
      ## make sure the dims consistant
      ## ----------------------------
      
      if(is.null(Y1) | is.null(Y2)){
        stop("Y1 or Y2 cannot be found\n")
      }
      if(any(is.na(Y1))){
        stop("NA values in Y1\n")
      }
      
      if(any(is.na(Y2))){
        stop("NA values in Y2\n")
      }
      Y1 = data.matrix(Y1)
      Y2 = data.matrix(Y2)
      
      # check dimensions of Y, Y1, Y2
      if(ncol(Y2) != nGene | nrow(Y2) != nSam |
         ncol(Y1) != nGene | nrow(Y1) != nSam)
        stop("dimensions of Y, Y1 and Y2 do not match\n")
    }else{
      Y2 = Y1 = Y
    }
    
    if(length(GeneSnpList) != 0 & length(GeneSnpList) !=  nGene){
      stop("dimensions of GeneSnpList does not correct\n")
    }
    
    RcppT_trecase_mtest(Y=Y, Y1=Y1, Y2=Y2, Z=Z, XX=XX, RHO=RHO, 
                        CNV1=CNV1, CNV2=CNV2, SNP_pos= SNPloc$pos, 
                        sChr=SNPloc$chr, gene_start=geneloc$start, 
                        gene_end=geneloc$end, gChr=geneloc$chr,
                        GeneSnpList=GeneSnpList, 
                        file_trec = file_trec, file_trecase = file_trecase, 
                        useLRT = useLRT, transTestP = transTestP, 
                        cis_window = cis_window, useASE = useASE, 
                        min_ASE_total = min_ASE_total, min_nASE = min_nASE, 
                        min_nASE_het = min_nASE_het, eps = eps, 
                        max_iter = max_iter, show = show) 
    
  }
