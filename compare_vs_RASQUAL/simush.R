args = commandArgs(trailingOnly = TRUE)
#args = c("0", "0", "0.5", "1", "500","2")
percase = 0.1
args
b0 = as.numeric(args[1])
th = as.numeric(args[2])
dv = as.numeric(args[3])
niter = as.numeric(args[4])
ss = as.numeric(args[5])
nsnp = as.numeric(args[6])
b1 = 0
#mn = 50
mn = 100
#rho = 0
set.seed(12345)
workdir = "/pine/scr/z/h/zhabotyn/R01/2019_03_20"
setwd(workdir)

phi = th/dv
rho = phi/(1+phi)
c(rho, th)

library(rcppreqtl, lib.loc="/nas/longleaf/home/zhabotyn/progs/RlibDev")
#library(VGAM)
#library(MASS)
ls("package:rcppreqtl")

Xm = makeXmatr(ss=2)
#betas = c(1, .8, .6, .4)#betas = c(1, .8, .6, .4, .2)
betas = c(3,.2,.05,.5)
#niter = 100; b0=0.5
if(nsnp==4){
  dat4 = simu4(niter, Xmatr=Xm$Xmatr, haplotype=Xm$thp, totmean=100, 
     percase = percase, dblcnt = 0, 
     phiNB = th, phiBB = th/dv, b0 = b0, b1 = 0, betas = betas) 
     dat0 = dat4; 
     dat0$asn = dat0$asnA; dat0$asnp =dat0$asnpA; 
     dat0$haplotype = dat0$haplotypeA
}else if(nsnp==2){
  dat2 = simu2(niter, Xmatr=Xm$Xmatr, haplotype=Xm$thp, totmean=100, 
     percase = percase, dblcnt = 0, 
     phiNB = th, phiBB = th/dv, b0 = b0, b1 = 0, betas = betas) 
     dat0 = dat2; 
     dat0$asn = dat0$asnA; dat0$asnp =dat0$asnpA; 
     dat0$haplotype = dat0$haplotypeA
}else if(nsnp==1){
   dat2 = simu2(niter, Xmatr=Xm$Xmatr, haplotype=Xm$thp, totmean=100, 
     percase = percase, dblcnt = 0, 
     phiNB = th, phiBB = th/dv, b0 = b0, b1 = 0, betas = betas) 
     dat0 = dat2; 
}

#dat0 = dat4
st3 = proc.time()
fm = fit1m(subset=1:niter, data=dat0)
en3 = proc.time()
en3[3]-st3[3]


suff = sprintf("add%s_phi%s_div%s_iter%s_ss%s_nsnp%s", 
               b0, th, dv, niter, ss, nsnp)
suff
if(!file.exists("resnew"))dir.create("resnew")

write.table(fm$full, sprintf("resnew/full_%s.csv", suff),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$short,sprintf("resnew/short_%s.csv", suff),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testadds,sprintf("resnew/short_testadd_%s.csv", suff),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testaddf,sprintf("resnew/full_testadd_%s.csv", suff),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(en3[3]-st3[3],sprintf("resnew/time_%s.csv", suff),
            quote=F, row.names=F, col.names=F, sep=",")

alpha = 0.05
c(mean(fm$full[,9]<alpha), mean(fm$short[,9]<alpha))



q("no")
