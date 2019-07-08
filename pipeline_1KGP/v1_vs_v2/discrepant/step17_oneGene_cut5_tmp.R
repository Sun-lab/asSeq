  tot = read.table("gene_tot.dat", as.is=T)
  a1 = read.table("gene_a1.dat", as.is=T)
  a2 = read.table("gene_a2.dat", as.is=T)
  inf = read.table("gene_inf.dat", as.is=T)
  gens = read.table("gene_gens.dat", as.is=T)
  vcfs = read.table("gene_vcfs.dat", as.is=T)
  sam = read.table("Xmat.dat", as.is=T, header=F)
  sam = matrix(unlist(sam), nrow=nrow(sam))
  
  tot = matrix(unlist(tot), nrow=1)
  a1 = matrix(unlist(a1), nrow=1)
  a2 = matrix(unlist(a2), nrow=1)
  
  nr=nrow(gens)
  gens = matrix(as.numeric(unlist(gens)),nrow=nr)
  chri = inf[1,2]
  #v1     
  library(asSeq, lib.loc="/nas/longleaf/home/zhabotyn/progs/Rlib")
  mChr = rep(chri, nr)
  eChr = rep(chri, nrow(tot))
  
  res.fil = sprintf("v1_%s", inf[,1])
  res.lon = sprintf("%s_eqtl.txt", res.fil)
  
  
  ePos = inf[,3]
  eEnd = inf[,4]
  eExt = eEnd-ePos
  pos = round((ePos+eEnd)/2)
  vcfi = vcfs[,1]
  posi = vcfs[,2]
  
  trecase(Y=t(tot), Y1=t(a1), Y2=t(a2), X=sam, Z=t(gens[,-1]), output.tag=res.fil, p.cut=1,
  eChr=eChr, ePos=pos, mChr=mChr, mPos=gens[,1], local.distance=2e5)
  
  eqtl = read.table(res.lon, header=T, as.is=T)
  eqtl[,1] = posi[eqtl[,2]]; colnames(eqtl)[1] = "Pos"
  eqtl[,2] = vcfi[eqtl[,2]]
  write.table(eqtl, res.lon, row.names=F, col.names=T, quote=F, sep="\t")
  message("done version 1")
  
  
  
  #v2
  res.fil = sprintf("v2_%s", inf[,1])
  res.trecase = sprintf("%s_trecase.txt", res.fil)
  res.trec = sprintf("%s_trec.txt", res.fil)
  
  for(geni in 2:ncol(gens)){
    kp = gens[,geni] == 3
    gens[kp, geni] = 2
    kp = gens[,geni] == 4
    gens[kp, geni] = 3
  }
  
  pos = round((ePos+eEnd)/2)
  geneloc = inf
  colnames(geneloc) = c("gene", "chr", "start", "end")
  SNPloc = cbind(sprintf("SNP%s",1:nrow(gens)),rep(inf[1,2],nrow(gens)), gens[,1])
  SNPloc = data.frame(SNPloc)
  colnames(SNPloc) = c("id", "")
  
  detach("package:asSeq")
  library(asSeq2, lib.loc="/nas/longleaf/home/zhabotyn/progs/Rlib/")
# (Y, Y1 = NULL, Y2 = NULL, Z, XX, SNPloc, geneloc, GeneSnpList = list(), 
#    fam_nb = T, file_trec = "trec.txt", file_trecase = "trecase.txt", 
#    cis_window = 100000L, useASE = 1L, min_ASE_total = 8L, min_nASE = 5L, 
#    min_nASE_het = 5L, eps = 5e-05, max_iter = 400L, show = FALSE) 
res = trecase(Y=t(tot), Y1=t(a1), Y2=t(a2), XX=sam, Z=t(gens[,-1]), min_ASE_total = 5, min_nASE = 5,
              file_trecase=res.trecase, file_trec=res.trec, geneloc=geneloc, SNPloc=SNPloc, cis_window=2e5)

eqtl2 = read.table(res.trecase, header=T, as.is=T)
eqtl2[,1] = posi[eqtl2[,2]]; colnames(eqtl2)[1] = "Pos"
eqtl2[,2] = vcfi[eqtl2[,2]]
write.table(eqtl2, res.trecase, row.names=F, col.names=T, quote=F, sep="\t")

teqtl2 = read.table(res.trec, header=T, as.is=T)
teqtl2[,1] = posi[teqtl2[,2]]; colnames(eqtl)[1] = "Pos"
teqtl2[,2] = vcfi[teqtl2[,2]]
write.table(teqtl2, res.trec, row.names=F, col.names=T, quote=F, sep="\t")

eqtl2
eqtl


q("no")
