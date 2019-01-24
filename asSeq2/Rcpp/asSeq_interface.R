Rcpp::sourceCpp('/fh/fast/sun_w/licai/_tumor_eQTL/GitHub/asSeq/asSeq2/Rcpp/asSeq_test1.cpp')

trecase <- function(Y, Y1, Y2, Z, X, SNPinfo, Geneinfo, fam_nb, file_out,
                    cis_window = 1e5, useASE = 1, 
                    min_ASE_total=8, min_nASE=10, eps=1e-5,
                    max_iter=4000L, show = F){
  # R wrapper function for TReCASE
  #   Args:
  #     Y:  A matrix of TReC values.
  #     Y1: A matrix of ASE values of allele 1.
  #     Y2: A matrix of ASE values of allele 2.
  #     Z:  A matrix of Genotype value taking 1 (AA), 2 (AB), 3 (BA), 4(BB).
  #     X:  Covariates in TReC regression.
  #     
  #     SNPinfo:  SNP name, chr, pos.
  #     Geneinfo: Gene name, chr, start, end.
  X = data.matrix(X)
  # Variable Check
  if(nrow(Y)!=nrow(Z)){
    stop("y and z have different lengths")
  }
  
  if(nrow(Y)!=nrow(Y1)||ncol(Y)!=ncol(Y1)){
    stop("y and y1 have different dimensions.\n")
  }
  
  if(nrow(Y1)!=nrow(Y2)||ncol(Y1)!=ncol(Y2)){
    stop("y1 and y2 have different dimensions.\n")
  }
  
  if(!all(Z %in% c(0,1,2,3))){
    stop("z must have values of 0, 1, 2, or 3 \n")
  }
  
  if(!all(Y>=0)){
    stop("Y must be non-negative")
  }
  
  if(!all(Y>=0)){
    stop("Y1 must be non-negative")
  }
  
  if(!all(Y2>=0)){
    stop("Y2 must be non-negative")
  }
  
  if(!all((Y1+Y2)<=Y)){
    stop("Y1+Y2 must be less than or equal to Y!")
  }
  

  
  Nchr = table(SNPinfo[,2])
  if(length(Nchr) > 1){
  # match chromosome and position 
  # gene_names = colnames(Y)
  # snp_names  = colnames()
  # if(!all(colnames((Y1) == gene_names) | !all(colnames(Y2) == gene_name))){
  #   stop("colnames of Y1, Y2 must match colnames Y!")
  # }
  # 
  # genematch = match(gene_names, Geneinfo[,1], nomatch=0L)
  # 
  # if(all(genematch == 0))
  #   stop("gene name does not match those in the gene infomation")
  # 
  # snpmatch = match(snp_names, SNPinfo[,1], nomatch=0L)
  # 
  # if(all(snpmatch == 0))
  #   stop("SNP name does not match those in the SNP infomation")
  # 
  }else{
    # geneInd    = which(gene[,2] == names(Nchr))
    # Geneinfo   = Geneinfo[geneInd, 2]
    Gene_start = Geneinfo[,3] - cis_window
    Gene_end   = Geneinfo[,4] + cis_window
    GStrec = list()
    for(gg in 1:nrow(Geneinfo)){
      # gg = 1        
      genenamei = as.character(Geneinfo[gg,1])
      cat(genenamei, date(), '\n')
      y  = Y[, gg]
      y1 = Y1[, gg]
      y2 = Y2[, gg]
      
      lgy1 = Rcpp_lgy_add_1(y)
      gene_snp_res = NULL
      for(ss in 1:nrow(SNPinfo)){
        # ss = 1
        if(SNPinfo[ss,3] > Gene_start[gg] & SNPinfo[ss,3] < Gene_end[gg]){
          zz =zz2 = Z[, ss]
          zz[which(zz == 2)] = 1
          zz[which(zz == 3)] = 2
          # table(zz)
          if(useASE){
            ni = y1 + y2
            ind = which(ni>=min_ASE_total)
            
            if(length(ind) < min_nASE){
              next
            }
              ni  = ni[ind]
              ni0 = y1
              ni0[zz2 == 3] = y2[zz2 == 3]
              ni0 = ni0[ind]
              
              zeta = zz2 %in% c(1,2)
              zeta = zeta[ind]
              lbc  = lchoose(ni, ni0)
              
              res_ase = Rcpp_ase(ni, ni0, zeta, lbc, max_iter, eps, show) 
            
              res_trecase = Rcpp_trecase(y, X, zz, fam_nb, lgy1, ni, ni0, zeta, 
                                         lbc, max_iter, eps, show) 

          }
          res_trec = Rcpp_trec(y, X, zz, fam_nb, lgy1, max_iter, eps, show)
          if(useASE)
            gene_snp_res = rbind(gene_snp_res, c(SNPname = SNPinfo[ss,1,drop =F],
                              trec = unlist(res_trec),
                             ase = unlist(res_ase), 
                             trecase = unlist(res_trecase)))
          else
            gene_snp_res =  rbind(gene_snp_res, c(SNPname = SNPinfo[ss,1,drop =F],
                                                  trec = unlist(res_trec)))
        }

      }
    GStrec[[genenamei]] = gene_snp_res
  }
  
  }
  return(GStrec)
}

trecase2 <- function(Y, Y1, Y2, Z, X, SNPinfo, Geneinfo, fam_nb, file_out,
                    cis_window = 1e5, useASE = 1, 
                    min_ASE_total=8, min_nASE=10, eps=1e-5,
                    max_iter=4000L, show = F){
  # R wrapper function for TReCASE
  #   Args:
  #     Y:  A matrix of TReC values.
  #     Y1: A matrix of ASE values of allele 1.
  #     Y2: A matrix of ASE values of allele 2.
  #     Z:  A matrix of Genotype value taking 1 (AA), 2 (AB), 3 (BA), 4(BB).
  #     X:  Covariates in TReC regression.
  #     
  #     SNPinfo:  SNP name, chr, pos.
  #     Geneinfo: Gene name, chr, start, end.
  X = data.matrix(X)
  # Variable Check
  if(nrow(Y)!=nrow(Z)){
    stop("y and z have different lengths")
  }
  
  if(nrow(Y)!=nrow(Y1)||ncol(Y)!=ncol(Y1)){
    stop("y and y1 have different dimensions.\n")
  }
  
  if(nrow(Y1)!=nrow(Y2)||ncol(Y1)!=ncol(Y2)){
    stop("y1 and y2 have different dimensions.\n")
  }
  
  if(!all(Z %in% c(0,1,2,3))){
    stop("z must have values of 0, 1, 2, or 3 \n")
  }
  
  if(!all(Y>=0)){
    stop("Y must be non-negative")
  }
  
  if(!all(Y>=0)){
    stop("Y1 must be non-negative")
  }
  
  if(!all(Y2>=0)){
    stop("Y2 must be non-negative")
  }
  
  if(!all((Y1+Y2)<=Y)){
    stop("Y1+Y2 must be less than or equal to Y!")
  }
  
  
  
  Nchr = table(SNPinfo[,2])
  if(length(Nchr) > 1){

    
  }else{
    # geneInd    = which(gene[,2] == names(Nchr))
    # Geneinfo   = Geneinfo[geneInd, 2]
    Gene_start = Geneinfo[,3] - cis_window
    Gene_end   = Geneinfo[,4] + cis_window
    GStrec = list()
    minPtrecase = NULL
    for(gg in 1:nrow(Geneinfo)){
      # gg = 1        
      ptmp = 1
      genenamei = as.character(Geneinfo[gg,1])
      cat(genenamei, date(), '\n')
      y  = Y[, gg]
      y1 = Y1[, gg]
      y2 = Y2[, gg]
      
      lgy1 = Rcpp_lgy_add_1(y)
      gene_snp_res = NULL
      for(ss in 1:nrow(SNPinfo)){
        # ss = 1
        if(SNPinfo[ss,3] > Gene_start[gg] & SNPinfo[ss,3] < Gene_end[gg]){
          zz =zz2 = Z[, ss]
          zz[which(zz == 2)] = 1
          zz[which(zz == 3)] = 2
          # table(zz)

          res_trec = Rcpp_trec(y, X, zz, fam_nb, lgy1, max_iter, eps, show)
          if(res_trec$pvalue < ptmp){
            ptmp = res_trec$pvalue
            res_trecMin = res_trec
            ssMin = ss
          }
            
          gene_snp_res =  rbind(gene_snp_res, c(SNPname = SNPinfo[ss,1,drop =F],
                                                  trec = unlist(res_trec)))
          
        }
      }
      
      zz =zz2 = Z[, ssMin]
      zz[which(zz == 2)] = 1
      zz[which(zz == 3)] = 2
      
      ni = y1 + y2
      ind = which(ni>=min_ASE_total)
      
      if(length(ind) < min_nASE){
        next
      }
      ni  = ni[ind]
      ni0 = y1
      ni0[zz2 == 3] = y2[zz2 == 3]
      ni0 = ni0[ind]
      
      zeta = zz2 %in% c(1,2)
      zeta = zeta[ind]
      lbc  = lchoose(ni, ni0)
      
      res_ase = Rcpp_ase(ni, ni0, zeta, lbc, max_iter, eps, show) 
      res_trecase = Rcpp_trecase(y, X, zz, fam_nb, lgy1, ni, ni0, zeta, 
                                 lbc, max_iter, eps, show) 
      
      minPtrecase =  rbind(minPtrecase, c(GeneName = Geneinfo[gg,1,drop=F],
                                    SNPname = SNPinfo[ssMin,1,drop =F],
                                    trec = unlist(res_trecMin),
                                    ase = unlist(res_ase), 
                                    trecase = unlist(res_trecase))) 
      
      GStrec[[genenamei]] = gene_snp_res
    }
     
  }
  return(list(GStrec, minPtrecase))
}

