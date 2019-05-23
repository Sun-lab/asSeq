#setwd("C:/Users/Vasyl/Documents/GitHub/asSeq/compare_vs_RASQUAL/fisher")
library(VGAM)
source("fishercl.R")


nas = 20
phi = .1
#phi = .5
for(phi in c(.1, .5)){
  rho = phi/(1+phi)
  par3 = par4 =
  par1 = par2 = matrix(0, nrow=nas, ncol=4)
  for(i in 1:nas){
  f8 = f4 = f2 = f1 = matrix(0,nrow=2, ncol=2)
  ni = i+(0);ni0=0
  for(j in 1:1){
    f8 = f8 + fishBB(n = rep(ni,8), y=rep(ni0,8), rho=rho, pi1=.5)
    f4 = f4 + fishBB(n = rep(2*ni, 4), y=rep(ni0,4), rho=rho, pi1=.5)
    f2 = f2 + fishBB(n = rep(4*ni, 2), y=rep(ni0,2), rho=rho, pi1=.5)
    f1 = f1 + fishBB(n = rep(8*ni, 1), y=rep(ni0,1), rho=rho, pi1=.5)
    ni = ni+1
    }
    par1[i,1] = f1[1,1]
    par1[i,2] = f2[1,1]
    par1[i,3] = f4[1,1]
    par1[i,4] = f8[1,1]
    
    par2[i,1] = f1[2,2]
    par2[i,2] = f2[2,2]
    par2[i,3] = f4[2,2]
    par2[i,4] = f8[2,2]
  
    for(j in 1:4){
    m = matrix(0, nrow=2, ncol=2)
    m[1,1] = par1[i,j]
    m[2,2] = par2[i,j]
    m1 = qr.solve(m)
    par3[i,j] = m1[1,1]
    par4[i,j] = m1[2,2]
    }
    message(i)
  }
  par3[1,4]=1e10
  if(phi==0.1){
    par1.1 = par1
    par2.1 = par2
    par3.1 = par3
    par4.1 = par4
  }else{
    par1.5 = par1
    par2.5 = par2
    par3.5 = par3
    par4.5 = par4
  }
}

cexes=1.3
cols = c("black", "blue", "goldenrod", "red")
phi = .1
png(sprintf("sd_by_asc_1ind_%s.png", phi), height=4, width=8, units="in", res=300)
par(mfrow=c(1,2))
ind1 = 1:nas
ylim = range(c(sqrt(par4.1[ind1,1]), sqrt(par4.1[ind1,2]), sqrt(par4.1[ind1,3]), sqrt(par4.1[ind1,4])))
plot(ind1, sqrt(par4.1[ind1,1]), xlim=c(1,nas), ylim=ylim, type="l", bty="n", bty="n", 
xlab="ASReC/SNP, 8 SNP case", ylab="sd", main="(a) eQTL(prop)", cex.main=cexes, cex.lab=cexes, cex.axis=cexes)
axis(1, at=c(1,5), labels=c(1,5), cex.lab=cexes, cex.axis=cexes)
abline(v=1, lty=3, col="grey")
for(i in 2:4)lines(ind1, sqrt(par4.1[ind1,i]), col=cols[i])
legend("topright", legend=c(sprintf("OD(theta)=%s",phi),sprintf("OD(rho)=%s",round(phi/(1+phi),3))), bty="n",cex=cexes) 

ylim = range(c(sqrt(par3.1[ind1,])),na.rm=T)
ylim[ylim>.5] = .5
plot(ind1, sqrt(par3.1[ind1,1]), xlim=c(1,nas), ylim=ylim, type="l", bty="n", 
xlab="ASReC/SNP, 8 SNP case", ylab="sd", main="(b) OD (rho)", cex.main=cexes, cex.lab=cexes, cex.axis=cexes)
axis(1, at=c(1,5), labels=c(1,5), cex.lab=cexes, cex.axis=cexes)
abline(v=1, lty=3, col="grey")
ylim = range(c(par3.1[ind1,1], par3.1[ind1,2], par3.1[ind1,3], par3.1[ind1,4]))
for(i in 2:4)lines(ind1, sqrt(par3.1[ind1,i]), col=cols[i])
legend("topright", legend=c("1 SNP", "2 SNP", "4 SNP", "8 SNP"), 
bty="n", text.col=cols,cex=cexes)
dev.off()

ind1 = 1:nas
phi = .5
png(sprintf("sd_by_asc_1ind_%s.png", phi), height=4, width=8, units="in", res=300)
par(mfrow=c(1,2))
ylim = range(c(sqrt(par4.5[ind1,1]), sqrt(par4.5[ind1,2]), sqrt(par4.5[ind1,3]), sqrt(par4.5[ind1,4])))
plot(ind1, sqrt(par4.5[ind1,1]), xlim=c(1,nas), ylim=ylim, type="l", bty="n", bty="n", 
xlab="ASReC/SNP, 8 SNP case", ylab="sd", main="(c) eQTL(prop)", cex.main=cexes, cex.lab=cexes, cex.axis=cexes)
axis(1, at=c(1,5), labels=c(1,5), cex.lab=cexes, cex.axis=cexes)
abline(v=1, lty=3, col="grey")
for(i in 2:4)lines(ind1, sqrt(par4.5[ind1,i]), col=cols[i])
legend("topright", legend=c(sprintf("OD(theta)=%s",phi),"",sprintf("OD(rho)=%s",round(phi/(1+phi),3))), bty="n",cex=cexes) 

ylim = range(c(sqrt(par3.5[ind1,])),na.rm=T)
ylim[ylim>.5] = .5
plot(ind1, sqrt(par3.5[ind1,1]), xlim=c(1,nas), ylim=ylim, type="l", bty="n", 
xlab="ASReC/SNP, 8 SNP case", ylab="sd", main="(d) OD (rho)", cex.main=cexes, cex.lab=cexes, cex.axis=cexes)
axis(1, at=c(1,5), labels=c(1,5), cex.lab=cexes, cex.axis=cexes)
abline(v=1, lty=3, col="grey")
ylim = range(c(par3.5[ind1,1], par3.5[ind1,2], par3.5[ind1,3], par3.5[ind1,4]))
for(i in 2:4)lines(ind1, sqrt(par3.5[ind1,i]), col=cols[i])
legend("topright", legend=c("1 SNP", "2 SNP", "4 SNP", "8 SNP"), 
bty="n", text.col=cols,cex=cexes)
dev.off()


phi = .1
png(sprintf("relSD_%s.png", phi), height=4, width=8, units="in", res=300)
par(mfrow=c(1,2))
y = par4.1/par4.1[,1]
ylim = c(0, 1)
plot(y[,1], xlim=c(1,nas), ylim=ylim, type="l", bty="n", bty="n", 
xlab="ASReC/SNP, 8 SNP case", ylab="relative sd", main="(a) eQTL(prop)", cex.main=cexes, cex.lab=cexes, cex.axis=cexes)
for(i in 2:4){lines(ind1, y[,i], col=cols[i])}
legend("topright", legend=c(sprintf("OD(theta)=%s",phi),sprintf("OD(rho)=%s",round(phi/(phi+1),3))), bty="n",cex=cexes)  
axis(1, at=c(1,5), labels=c(1,5), cex.lab=cexes, cex.axis=cexes)

y = par3.1/par3.1[,1]
ylim = c(0,3.5)
plot(y[,1], xlim=c(1,nas), ylim=ylim, type="l", bty="n", bty="n", 
xlab="ASReC/SNP, 8 SNP case", ylab="relative sd", main="(b) OD(rho)", cex.main=cexes, cex.lab=cexes, cex.axis=cexes)
for(i in 2:4){lines(ind1, y[,i], col=cols[i])}
legend("topright", legend=c("1 SNP", "2 SNP", "4 SNP", "8 SNP"), 
bty="n", text.col=cols,cex=cexes) 
axis(1, at=c(1,5), labels=c(1,5), cex.lab=cexes, cex.axis=cexes)

dev.off()

phi = .5
png(sprintf("relSD_%s.png", phi), height=4, width=8, units="in", res=300)
par(mfrow=c(1,2))
y = par4.5/par4.5[,1]
ylim = c(0, 1)
plot(y[,1], xlim=c(1,nas), ylim=ylim, type="l", bty="n", bty="n", 
xlab="ASReC/SNP, 8 SNP case", ylab="relative sd", main="(c) eQTL(prop)", cex.main=cexes, cex.lab=cexes, cex.axis=cexes)
for(i in 2:4){lines(ind1, y[,i], col=cols[i])}
legend("topright", legend=c(sprintf("OD(theta)=%s",phi),sprintf("OD(rho)=%s",round(phi/(phi+1),3))), bty="n",cex=cexes)  
axis(1, at=c(1,5), labels=c(1,5), cex.lab=cexes, cex.axis=cexes)

y = par3.5/par3.5[,1]
ylim = c(0,3.5)
plot(y[,1], xlim=c(1,nas), ylim=ylim, type="l", bty="n", bty="n", 
xlab="ASReC/SNP, 8 SNP case", ylab="relative sd", main="(d) OD(rho)", cex.main=cexes, cex.lab=cexes, cex.axis=cexes)
for(i in 2:4){lines(ind1, y[,i], col=cols[i])}
legend("topright", legend=c("1 SNP", "2 SNP", "4 SNP", "8 SNP"), 
bty="n", text.col=cols,cex=cexes) 
axis(1, at=c(1,5), labels=c(1,5), cex.lab=cexes, cex.axis=cexes)

dev.off()

png(sprintf("hess_log10scale_%s.png", phi), height=4, width=8, units="in", res=300)
par(mfrow=c(1,2))
x1 = log10(par1[ind1,1]); x1[x1 < 0] = 0
y1 = log10(par1[ind1,4]); y1[y1 < 0] = 0;y1[!is.finite(y1)]=0
x2 = log10(par1[ind1,1]); x2[x2 < 0] = 0
y2 = log10(par1[ind1,3]); y2[y2 < 0] = 0
x3 = log10(par1[ind1,1]); x3[x3 < 0] = 0
y3 = log10(par1[ind1,2]); y3[y3 < 0] = 0
xlim=range(c(x1,x2,x3))
ylim=range(c(y1,y2,y3))
plot(x1, y1, type="l", col="red", bty="n", main="OD(rho)", xlab="1 SNP", ylab="multi SNP", xlim=xlim,ylim=ylim) 
lines(x2, y2, col="goldenrod") 
lines(x3, y3, col="blue") 
abline(a=0, b=1, lty=3)
legend("topleft", legend=c("2 SNP", "4 SNP", "8 SNP"), 
bty="n", text.col=c("blue", "goldenrod", "red"))

x1 = log10(par2[ind1,1]); x1[x1 < 0] = 0
y1 = log10(par2[ind1,4]); y1[y1 < 0] = 0;y1[!is.finite(y1)]=0
x2 = log10(par2[ind1,1]); x2[x2 < 0] = 0
y2 = log10(par2[ind1,3]); y2[y2 < 0] = 0
x3 = log10(par2[ind1,1]); x3[x3 < 0] = 0
y3 = log10(par2[ind1,2]); y3[y3 < 0] = 0
xlim=range(c(x1,x2,x3))
ylim=range(c(y1,y2,y3))
plot(x1, y1, type="l", col="red", ylim=ylim, bty="n", main="eQTL(prop)", xlab="1 SNP", ylab="multi SNP") 
lines(x2, y2, col="goldenrod") 
lines(x3, y3, col="blue") 
abline(a=0, b=1, lty=3)
legend("topleft", legend=c("2 SNP", "4 SNP", "8 SNP"), 
bty="n", text.col=c("blue", "goldenrod", "red"))
dev.off()

png(sprintf("hess_relation_%s.png", phi), height=4, width=8, units="in", res=300)
par(mfrow=c(1,2))
x1 = par1[ind1,1]; x1[x1 < 0] = 0
y1 = par1[ind1,4]; y1[y1 < 0] = 0
x2 = par1[ind1,1]; x2[x2 < 0] = 0
y2 = par1[ind1,3]; y2[y2 < 0] = 0
x3 = par1[ind1,1]; x3[x3 < 0] = 0
y3 = par1[ind1,2]; y3[y3 < 0] = 0
ylim=range(c(y1,y2,y3))
plot(x1, y1, type="l", col="red", main="OD(rho)", xlab="1 SNP", ylab="multi SNP", bty="n", ylim=ylim) 
lines(x2, y2, col="goldenrod") 
lines(x3, y3, col="blue") 
abline(a=0, b=1, lty=3)
legend("topleft", legend=c("2 SNP", "4 SNP", "8 SNP"), 
bty="n", text.col=c("blue", "goldenrod", "red"))

x1 = par2[ind1,1]; x1[x1 < 0] = 0
y1 = par2[ind1,4]; y1[y1 < 0] = 0
x2 = par2[ind1,1]; x2[x2 < 0] = 0
y2 = par2[ind1,3]; y2[y2 < 0] = 0
x3 = par2[ind1,1]; x3[x3 < 0] = 0
y3 = par2[ind1,2]; y3[y3 < 0] = 0
ylim=range(c(y1,y2,y3))
plot(x1, y1, type="l", col="red", bty="n", main="eQTL(prop)", xlab="1 SNP", ylab="multi SNP", ylim=ylim)
lines(x2, y2, col="goldenrod") 
lines(x3, y3, col="blue") 
abline(a=0, b=1, lty=3)
legend("topleft", legend=c("2 SNP", "4 SNP", "8 SNP"), 
bty="n", text.col=c("blue", "goldenrod", "red"))
dev.off()


cols=c("black", "blue", "goldenrod", "red")
png(sprintf("hess_logscale1_%s.png",phi), height=4, width=8, units="in", res=300)
par(mfrow=c(1,2))
ylim=range(c(log10(par1[ind1,])), na.rm=T);ylim[ylim<0]=0
i = 1
x = log10(par1[ind1,i]); x[x < 0] = 0
plot(ind1, x, type="l", col=cols[i], bty="n", main="OD(rho)", ylab="information", xlab="ASReC/SNP @ 8 SNP scenario", ylim=ylim, xlim=c(1,20)) 
axis(1, at=c(1,5), labels=c(1,5))
abline(v=1, lty=3, col="grey")
for(i in 2:4){
x = log10(par1[ind1,i]); x[x < 0] = 0
lines(ind1, x, col=cols[i]) 
}
legend("bottomright", legend=c("1 SNP", "2 SNP", "4 SNP", "8 SNP"), bty="n", text.col=cols)

i = 1
ylim=range(c(log10(par2[ind1,])), na.rm=T)
x = log10(par2[ind1,i]); x[x < 0] = 0
plot(ind1, x, type="l", col=cols[i], bty="n", main="eQTL(prop)", ylab="information", xlab="ASReC/SNP @ 8 SNP scenario", ylim=ylim, xlim=c(1,20)) 
axis(1, at=c(1,5), labels=c(1,5))
abline(v=1, lty=3, col="grey")
for(i in 2:4){
x = log10(par2[ind1,i]); x[x < 0] = 0
lines(ind1, x, col=cols[i]) 
}
legend("topleft", sprintf("OD(theta)=%s", phi), bty="n")
legend("bottomright", legend=c("1 SNP", "2 SNP", "4 SNP", "8 SNP"), bty="n", text.col=cols)
dev.off()

png(sprintf("hess_relation1_%s.png", phi), height=4, width=8, units="in", res=300)
par(mfrow=c(1,2))
i = 1
ylim=range(c((par1[ind1,])), na.rm=T)
x = (par1[ind1,i]); x[x < 0] = 0
plot(ind1, x, type="l", col=cols[i], bty="n", main="OD(rho)", ylab="information", xlab="ASReC/SNP @ 8 SNP scenario", ylim=ylim, xlim=c(1,20)) 
axis(1, at=c(1,5), labels=c(1,5))
abline(v=1, lty=3, col="grey")
for(i in 2:4){
x = (par1[ind1,i]); x[x < 0] = 0
lines(ind1, x, col=cols[i]) 
}
legend("topleft", legend=c("1 SNP", "2 SNP", "4 SNP", "8 SNP"), bty="n", text.col=cols)

i = 1
ylim=range(c((par2[ind1,])), na.rm=T)
x = (par2[ind1,i]); x[x < 0] = 0
plot(ind1, x, type="l", col=cols[i], bty="n", main="eQTL(prop)", ylab="information", xlab="ASReC/SNP @ 8 SNP scenario", ylim=ylim, xlim=c(1,20)) 
axis(1, at=c(1,5), labels=c(1,5))
abline(v=1, lty=3, col="grey")
for(i in 2:4){
x = (par2[ind1,i]); x[x < 0] = 0
lines(ind1, x, col=cols[i]) 
}
legend("topleft", sprintf("OD(theta)=%s", phi), bty="n")
#legend("bottomright", legend=c("1 SNP", "2 SNP", "4 SNP", "8 SNP"), bty="n", text.col=cols)
dev.off()


q("no")
