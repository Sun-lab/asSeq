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
  #detach("package:asSeq2")
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
  eChr=eChr, ePos=pos, mChr=mChr, mPos=gens[,1], local.distance=2e5, trace=2)
  
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

SNPloc = data.frame(snp=sprintf("SNP%s",1:nrow(gens)),
                    chr=rep(inf[1,2],nrow(gens)),
                    pos=gens[,1], stringsAsFactors = F)
str(SNPloc)
  
  detach("package:asSeq")
  library(asSeq2, lib.loc="/nas/longleaf/home/zhabotyn/progs/Rlib/")
res = trecase(Y=t(tot), Y1=t(a1), Y2=t(a2), XX=sam, Z=t(gens[,-1]), min_ASE_total = 5, min_nASE = 5,
              file_trecase=res.trecase, file_trec=res.trec, geneloc=geneloc, SNPloc=SNPloc, cis_window=2e5)


eqtl2 = read.table(res.trecase, header=T, as.is=T)
eqtl2[,1] = posi[eqtl2[,2]]; colnames(eqtl2)[1] = "Pos"
eqtl2[,2] = vcfi[eqtl2[,2]]
write.table(eqtl2, res.trecase, row.names=F, col.names=T, quote=F, sep="\t")


eqtl2[1:4,]
eqtl[1:4,]

log10(eqtl2$Joint_Pvalue)-log10(eqtl$Joint_Pvalue)

q("no")
