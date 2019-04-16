#module add tabix;module add gcc;module add gsl
args = commandArgs(trailingOnly = TRUE)
#args = c("0", "0.1", "1", "10", "2")
#args = c("1", "0.1", "1", "10", "2")
percase = 0.1
args
b0 = as.numeric(args[1])
th = as.numeric(args[2])
dv = as.numeric(args[3])
niter = as.numeric(args[4])
ss = as.numeric(args[5])
#mn = 50
mn = 100
nsnp = 8
set.seed(12345)

phi = th/dv
rho = phi/(1+phi)
c(rho, th)

  library(rcppreqtl)
  ls("package:rcppreqtl")
  
  Xm = makeXmatr(ss=ss)
  betas = c(3,.2,.05,.5)

  #simu8 the data according to provided covariates
  #in addition to one allele-specific count per individual it produces
  #the allele-specific reads equally split among 8, 4 and 2 SNPS
  
    dat8 = simu8(niter, Xmatr=Xm$Xmatr, haplotype=Xm$thp, totmean=mn, dblcnt = 0, 
       percase = percase, phiNB = th, phiBB = th/dv, b0 = b0, b1 = 0, betas = betas)

  nind = ncol(dat8$haplotype)
  hap = c("0|0", "0|1", "1|0", "1|1")

  hap1 = "0|1"
  hap2 = "1|0"
  
  hapr = matrix(sprintf("%s:%s,%s", hap[dat8$haplotype+1], 0, 0), nrow=niter)
  ind1 = 1:nind
  ind2 = 1:nind+nind
  ind3 = 1:nind+nind*2
  ind4 = 1:nind+nind*3
  ind5 = 1:nind+nind*4
  ind6 = 1:nind+nind*5
  ind7 = 1:nind+nind*6
  ind8 = 1:nind+nind*7

  #rewrite for RASQUAL
  #rewriting in VCF format for RASQUAL - we switch reference and alternative among SNPs    
  ind=ind1; hap1S = matrix(sprintf("%s:%s,%s", hap1, 
                 (dat8$asn8-dat8$asnp8)[,ind], dat8$asnp8[,ind]), nrow=niter)
  ind=ind2; hap2S = matrix(sprintf("%s:%s,%s", hap2, 
                 dat8$asnp8[,ind],(dat8$asn8-dat8$asnp8)[,ind]),  nrow=niter)
  ind=ind3; hap3S = matrix(sprintf("%s:%s,%s", hap1, 
                 (dat8$asn8-dat8$asnp8)[,ind], dat8$asnp8[,ind]), nrow=niter)
  ind=ind4; hap4S = matrix(sprintf("%s:%s,%s", hap2, 
                 dat8$asnp8[,ind],(dat8$asn8-dat8$asnp8)[,ind]),  nrow=niter)
  ind=ind5; hap5S = matrix(sprintf("%s:%s,%s", hap1, 
                 (dat8$asn8-dat8$asnp8)[,ind], dat8$asnp8[,ind]), nrow=niter)
  ind=ind6; hap6S = matrix(sprintf("%s:%s,%s", hap2, 
                 dat8$asnp8[,ind],(dat8$asn8-dat8$asnp8)[,ind]),  nrow=niter)
  ind=ind7; hap7S = matrix(sprintf("%s:%s,%s", hap1, 
                 (dat8$asn8-dat8$asnp8)[,ind], dat8$asnp8[,ind]), nrow=niter)
  ind=ind8; hap8S = matrix(sprintf("%s:%s,%s", hap2, 
                 dat8$asnp8[,ind],(dat8$asn8-dat8$asnp8)[,ind]),  nrow=niter)

  #adding corresponding VCF info    
  chr = 1
  pos = posr = (1:niter)*1e1
  pos1 = (1:niter)*1e1+1
  pos2 = (1:niter)*1e1+2
  pos3 = (1:niter)*1e1+3
  pos4 = (1:niter)*1e1+4
  pos5 = (1:niter)*1e1+5
  pos6 = (1:niter)*1e1+6
  pos7 = (1:niter)*1e1+7
  pos8 = (1:niter)*1e1+8
  rs = rsr = sprintf("rs%s", 1:niter)
  
  rs1 = sprintf("%sf1", rs)
  rs2 = sprintf("%sf2", rs)
  rs3 = sprintf("%sf3", rs)
  rs4 = sprintf("%sf4", rs)
  rs5 = sprintf("%sf5", rs)
  rs6 = sprintf("%sf6", rs)
  rs7 = sprintf("%sf7", rs)
  rs8 = sprintf("%sf8", rs)
  
  bs1 = rep("A", niter)
  bs2 = rep("C", niter)
  oth = rep(100, niter)
  
  fil = rep("PASS", niter)
  rsq = sprintf("RSQ=%s",rep(1, niter))
  inf = rep("GT:AS", niter)

  hapmr = data.frame(chr, pos=posr, rs=rsr, bs1, bs2, oth, fil, rsq, inf, hapr)
  hapmf1 = data.frame(chr, pos1, rs1, bs1, bs2, oth, fil, rsq, inf, hap1S)
  hapmf2 = data.frame(chr, pos2, rs2, bs2, bs1, oth, fil, rsq, inf, hap2S)
  hapmf3 = data.frame(chr, pos3, rs3, bs1, bs2, oth, fil, rsq, inf, hap3S)
  hapmf4 = data.frame(chr, pos4, rs4, bs2, bs1, oth, fil, rsq, inf, hap4S)
  hapmf5 = data.frame(chr, pos5, rs5, bs1, bs2, oth, fil, rsq, inf, hap5S)
  hapmf6 = data.frame(chr, pos6, rs6, bs2, bs1, oth, fil, rsq, inf, hap6S)
  hapmf7 = data.frame(chr, pos7, rs7, bs1, bs2, oth, fil, rsq, inf, hap7S)
  hapmf8 = data.frame(chr, pos8, rs8, bs2, bs1, oth, fil, rsq, inf, hap8S)
  colnames(hapmf1) = colnames(hapmr)                                    
  colnames(hapmf2) = colnames(hapmr)                                    
  colnames(hapmf3) = colnames(hapmr)
  colnames(hapmf4) = colnames(hapmr)
  colnames(hapmf5) = colnames(hapmr)
  colnames(hapmf6) = colnames(hapmr)
  colnames(hapmf7) = colnames(hapmr)
  colnames(hapmf8) = colnames(hapmr)
  hapmf = rbind(hapmr, hapmf1, hapmf2, hapmf3, hapmf4, hapmf5, hapmf6, hapmf7, hapmf8)
  hapmf = hapmf[order(hapmf[,"pos"]),]
  
  #write each component RASQUAL style - VCF, offset matrix K, covariate matrix X and total counts - binary files
  options("scipen"=100, "digits"=4)
  if(!file.exists("rasq"))dir.create("rasq")
  suff = sprintf("add%s_phi%s_div%s_iter%s_ss%s_nsnp%s", b0, th, dv, niter, ss, nsnp)
  vcfn = sprintf("rasq/hap_%s.vcf", suff)
  totn = sprintf("rasq/tot_%s.txt", suff)
  totb = sprintf("rasq/tot_%s.bin", suff)
  vcfgz = sprintf("%s.gz", vcfn)
  kbin = sprintf("rasq/K_%s", suff)
  xbin = sprintf("rasq/X_%s", suff)
  
  write.table(hapmf,vcfn,quote=F,row.names=F,col.names=F, sep="\t")
  if(file.exists(vcfgz))system(sprintf("rm %s", vcfgz))
  system(sprintf("bgzip %s", vcfn))
  system(sprintf("tabix -p vcf %s", vcfgz))
  write.table(dat8$trc, totn, quote=F, col.names=F, sep="\t")
  writeBin(c(t(dat8$trc)), con=totb, double())
  writeBin(c((dat8$X)), con=xbin, double())
  
  kbinm =  exp(dat8$X[,4])
  kbinm = kbinm/median(kbinm)
  kbinm = matrix(rep(kbinm, niter), nrow=niter, byrow=T)
  writeBin(c(t(kbinm)), con=kbin, double())
  
  snps = sprintf("%s:%s-%s", chr, pos-1, pos+8)
  
  lrt0 = lrt = pis = rep(NA, niter)
  pim = cbind(pis, pis, pis)
  i = 1
  nms = c("fID","sID","SNPchr","SNPpos","Ref","Alt","AF","HWEChisq","impQual","qval",
  "Chi2","Pi","DeltaErr","PhiBias","OD","SNPID","NfSNP","NtestSNP","niterNull","niterAlt","rndtieloc","llh0","conv0","corfSNP", "corrSNP")
  
  tmpi = sprintf("rasq/tmp_%s.txt", suff)
  outi = sprintf("rasq/fin_%s.csv", suff)
  outitot = sprintf("rasq/fin_%s_tot.csv", suff)
  outiase = sprintf("rasq/fin_%s_ase.csv", suff)

  #start simulation gene by gene
  #in this simulation we consider fixed delta and phi and use all SNPs
  # 
  write.table(t(nms), outi, append=F, sep=",", quote=F, row.names=F, col.names=F)
  write.table(t(nms), outitot, append=F, sep=",", quote=F, row.names=F, col.names=F)
  write.table(t(nms), outiase, append=F, sep=",", quote=F, row.names=F, col.names=F)
  
  i = 0
  #niter = 10
  for(i in 1:niter){
  #for(i in 1:20){
  #i = i + 1
    snps = sprintf("%s:%s-%s", chr, max(1,pos[i]-1), pos[i]+9)

    rasBASE = sprintf("rasqual -y %s -k %s -x %s -p 4 -n %s", totb, kbin, xbin, nind)
    rasPOSj = sprintf("-j %s -l %s -m %s -s %s -e %s",i, max(pos8), max(pos8), 6, max(pos8))
    rasOPT = "-z -t --fix-delta --fix-phi --min-coverage-depth 0.0"
    cmdi = sprintf("tabix %s %s | %s %s -f %s %s --rsnp %s > %s",
      vcfgz, snps, rasBASE, rasPOSj, rownames(dat8$trc)[i], rasOPT, rsr[i], tmpi)
    system(cmdi)  
    gni = read.table(tmpi, header=F)
    write.table(gni, outi, append=T, sep=",", quote=F, row.names=F, col.names=F)
    
    rasOPTt = "-z -t --fix-delta --fix-phi --min-coverage-depth 0.0 --population-only"
    rasPOSt = sprintf("-j %s -l %s -m %s -s %s -e %s",i, max(pos8), max(pos8), pos[i+1]-1, pos[i+1]+1)
    cmdi = sprintf("tabix %s %s | %s %s -f %s %s --rsnp %s  > %s", 
      vcfgz, snps, rasBASE, rasPOSt, rownames(dat8$trc)[i], rasOPTt, rsr[i], tmpi)
    system(cmdi)  
    gni = read.table(tmpi, header=F)
    write.table(gni, outitot, append=T, sep=",", quote=F, row.names=F, col.names=F)
  
    rasOPTa = "-z -t --fix-delta --fix-phi --min-coverage-depth 0.0 --population-only"
    cmdi = sprintf("tabix %s %s | %s %s -f %s %s --rsnp %s  > %s",   
      vcfgz, snps, rasBASE, rasPOSj, rownames(dat8$trc)[i], rasOPTa, rsr[i], tmpi)
    system(cmdi)  
    gni = read.table(tmpi, header=F)
    write.table(gni, outiase, append=T, sep=",", quote=F, row.names=F, col.names=F)
  
    if(i%%10==0)message(i,"/",niter)
  }

  #
  file.remove(tmpi)
  
  
  
  jtfit = read.csv(outi)
  trcfit = read.csv(outitot)
  asefit = read.csv(outiase)
  
  jtpi = jtfit$Pi
  trcpi = trcfit$Pi
  asepi = asefit$Pi
  jtpi[jtfit$Ref=="C"] = 1- jtpi[jtfit$Ref=="C"]
  trcpi[trcfit$Ref=="C"] = 1- trcpi[trcfit$Ref=="C"]
  asepi[asefit$Ref=="C"] = 1- asepi[asefit$Ref=="C"]
  summary(jtpi)
  summary(trcpi)
  summary(asepi)
  1/(1+exp(-b0))                                
  
  summary(jtfit$OD)
  summary(trcfit$OD)
  summary(asefit$OD)
  
  
  
  summary(jtfit$Chi2)
  summary(trcfit$Chi2)
  summary(asefit$Chi2)


q("no")













q("no")
