args = commandArgs(trailingOnly = TRUE)
#args = c("0", "0.5", "1", "500","2", 4)
percase = 0.1
args
b0 = as.numeric(args[1])
th = as.numeric(args[2])
dv = as.numeric(args[3])
niter = as.numeric(args[4])
ss = as.numeric(args[5])
nsnp = as.numeric(args[6])
b1 = 0
mn = 100
set.seed(12345)

phi = th/dv
rho = phi/(1+phi)
c(rho, th)

library(rcppreqtl)
ls("package:rcppreqtl")

Xm = makeXmatr(ss=2)
betas = c(3,.2,.05,.5)

#simu4 and simu2 simulate the data according to provided covariates
#in addition to one allele-specific count per individual
#simu2 produces allele-specific counts split in 2 SNPs
#simu4 produces allele-specific split in 4 SNPs
#those additional fields are stored in asnA, asnpA and haplotypeA
if(nsnp==4){
  dat4 = simu4(niter, Xmatr=Xm$Xmatr, haplotype=Xm$thp, totmean=100, dblcnt = 0, 
     percase = percase, phiNB = th, phiBB = th/dv, b0 = b0, b1 = 0, betas = betas)
     #dat0 - data split into 4 SNPS in each individual
     #dat1 - one count per individual
     dat0 = dat4; 
     dat0$asn = dat0$asnA; dat0$asnp =dat0$asnpA; 
     dat0$haplotype = dat0$haplotypeA
     dat1 = dat4;
}else if(nsnp==2){
  dat2 = simu2(niter, Xmatr=Xm$Xmatr, haplotype=Xm$thp, totmean=100, dblcnt = 0,  
     percase = percase, phiNB = th, phiBB = th/dv, b0 = b0, b1 = 0, betas = betas)  
     #dat0 - data split into 4 SNPS in each individual
     #dat01 - one count per individual     
     dat0 = dat2; 
     dat0$asn = dat0$asnA; dat0$asnp =dat0$asnpA; 
     dat0$haplotype = dat0$haplotypeA
     dat1 = dat2;
}


alpha = 0.05
suff = sprintf("add%s_phi%s_div%s_iter%s_ss%s_nsnp%s", 
               b0, th, dv, niter, ss, nsnp)
suff

#fitting the data in multi-SNP scenario
st3 = proc.time()
fm = fit1m(subset=1:niter, data=dat0)
en3 = proc.time()
en3[3]-st3[3]
  
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
message(mean(fm$full[,9]<alpha), " ", mean(fm$short[,9]<alpha))





#fitting the data in one-SNP scenario
st3 = proc.time()
fm = fit1m(subset=1:niter, data=dat1)
en3 = proc.time()
en3[3]-st3[3]

if(!file.exists("resnew1"))dir.create("resnew1")
write.table(fm$full, sprintf("resnew1/full_%s.csv", suff),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$short,sprintf("resnew1/short_%s.csv", suff),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testadds,sprintf("resnew1/short_testadd_%s.csv", suff),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testaddf,sprintf("resnew1/full_testadd_%s.csv", suff),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(en3[3]-st3[3],sprintf("resnew1/time_%s.csv", suff),
            quote=F, row.names=F, col.names=F, sep=",")

message(mean(fm$full[,9]<alpha), " ", mean(fm$short[,9]<alpha))

q("no")
