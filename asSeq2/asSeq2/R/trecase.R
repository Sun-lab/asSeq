trecase <-
  function(Y, Y1, Y2, Z, XX, SNPloc, geneloc, fam_nb = T,
           file_trec = "trec.txt", file_trecase = "trecase.txt",
           cis_widow = 1e5L, useASE = 1L, min_ASE_total = 8L, min_nASE = 10L,
           eps = 5e-5, max_iter = 400L, show = FALSE)
  {

    minVar = 1e-8

    if(min_nASE < min_ASE_total){
      stop("min_nASE should not be smaller than min_ASE_total\n")
    }

    ##
    # check the NAs in X, Y, and Z
    #

    if(any(is.na(Y))){
      stop("NA values in Y\n")
    }

    if(any(is.na(XX))){
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

    Y = data.matrix(Y)
    Y1 = data.matrix(Y1)
    Y2 = data.matrix(Y2)
    Z = data.matrix(Z)
    XX = data.matrix(XX)


    nGene = ncol(Y)
    nSam  = nrow(Y)

    # check dimensions of Y, Y1, Y2
    if(ncol(Y2) != nGene | nrow(Y2) != nSam |
       ncol(Y1) != nGene | nrow(Y1) != nSam)
      stop("dimensions of Y, Y1 and Y2 do not match\n")

    # check dimension of Z
    if(nrow(Z) != nSam)
      stop("Z and Y have different sample size")

    # check X, Y, Z

    varY = apply(Y, 2, var)
    wVar = which(varY < minVar)

    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Y have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))

      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
      }


    varX  = apply(XX, 2, var)
    wVarX = which(varX < minVar)
    if(any(wVarX > 1)){
      stop("X has tiny variance")
    }

    varZ = apply(Z, 2, var)
    wVar = which(varZ < minVar)

    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Z have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))

      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
      }

    if(! all(as.numeric(Z) %in% c(0,1,2,3))){
      stop("Z must take values 0, 1, 2, or 3\n")
    }

    ##
    # check the length of chromosome location information
    #

    if(nrow(geneloc) != nGene)
      stop("number of Gene in GeneInfo doesn't match that in Y")
    if(nrow(SNPloc) != ncol(Z))
      stop("number of Gene in GeneInfo doesn't match that in Y")

    # Zh = Z
    # Z[Zh==2] = 1
    # Z[Zh==3] = 2

    Rcpp_trecase_mtest(Y, Y1, Y2, Z, XX, as.numeric(SNPloc[,3]),
              as.numeric(SNPloc[,2]), fam_nb, as.numeric(geneloc[,3]),
              as.numeric(geneloc[,4]), as.numeric(geneloc[,2]),
              file_trec, file_trecase,
              cis_widow, useASE, min_ASE_total, min_nASE,
              eps, max_iter, show)

  }

