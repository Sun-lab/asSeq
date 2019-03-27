args = commandArgs(trailingOnly = TRUE)
#args = c("0", "0.5", "1", "5","2", 4)
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

#simu8 the data according to provided covariates
#in addition to one allele-specific count per individual it produces
#the allele-specific reads equally split among 8, 4 and 2 SNPS

  dat8 = simu8(niter, Xmatr=Xm$Xmatr, haplotype=Xm$thp, totmean=100, dblcnt = 0, 
     percase = percase, phiNB = th, phiBB = th/dv, b0 = b0, b1 = 0, betas = betas)

     #dat0 - data split into 4 SNPS in each individual
     #dat1 - one count per individual

alpha = 0.05
if(!file.exists("resnew"))dir.create("resnew")
suff = sprintf("add%s_phi%s_div%s_iter%s_ss%s_nsnp", 
               b0, th, dv, niter, ss)
suff

#fitting the data in 8-SNP scenario
dat0 = dat8
dat0$asn = dat0$asn8S
dat0$asnp = dat0$asnp8S
dat0$haplotype = dat0$haplotype8S
nSNP = 8
st3 = proc.time()
fm = fit1m(subset=1:niter, data=dat0)
en3 = proc.time()
en3[3]-st3[3]
  
write.table(fm$full, sprintf("resnew/full_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$short,sprintf("resnew/short_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testadds,sprintf("resnew/short_testadd_%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testaddf,sprintf("resnew/full_testadd_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(en3[3]-st3[3],sprintf("resnew/time_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
message(mean(fm$full[,9]<alpha), " ", mean(fm$short[,9]<alpha))


#fitting the data in 4-SNP scenario
dat0 = dat8
dat0$asn = dat0$asn4S
dat0$asnp = dat0$asnp4S
dat0$haplotype = dat0$haplotype4S
nSNP = 4
st3 = proc.time()
fm = fit1m(subset=1:niter, data=dat0)
en3 = proc.time()
en3[3]-st3[3]
  
write.table(fm$full, sprintf("resnew/full_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$short,sprintf("resnew/short_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testadds,sprintf("resnew/short_testadd_%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testaddf,sprintf("resnew/full_testadd_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(en3[3]-st3[3],sprintf("resnew/time_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
message(mean(fm$full[,9]<alpha), " ", mean(fm$short[,9]<alpha))

#fitting the data in 2-SNP scenario
dat0 = dat8
dat0$asn = dat0$asn2S
dat0$asnp = dat0$asnp2S
dat0$haplotype = dat0$haplotype2S
nSNP = 2
st3 = proc.time()
fm = fit1m(subset=1:niter, data=dat0)
en3 = proc.time()
en3[3]-st3[3]
  
write.table(fm$full, sprintf("resnew/full_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$short,sprintf("resnew/short_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testadds,sprintf("resnew/short_testadd_%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testaddf,sprintf("resnew/full_testadd_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(en3[3]-st3[3],sprintf("resnew/time_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
message(mean(fm$full[,9]<alpha), " ", mean(fm$short[,9]<alpha))



#fitting the data in 1-SNP scenario
dat0 = dat8
nSNP = 1
st3 = proc.time()
fm = fit1m(subset=1:niter, data=dat0)
en3 = proc.time()
en3[3]-st3[3]
  
write.table(fm$full, sprintf("resnew/full_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$short,sprintf("resnew/short_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testadds,sprintf("resnew/short_testadd_%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(fm$testaddf,sprintf("resnew/full_testadd_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
write.table(en3[3]-st3[3],sprintf("resnew/time_%s%s.csv", suff, nSNP),
            quote=F, row.names=F, col.names=F, sep=",")
message(mean(fm$full[,9]<alpha), " ", mean(fm$short[,9]<alpha))




q("no")
