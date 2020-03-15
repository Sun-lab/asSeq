#dirs = c("run_Brain_Caudate_basal_ganglia_194_2e+05_long",
#"run_Brain_Caudate_basal_ganglia_194_2e+05_short",
#"run_Brain_Caudate_basal_ganglia_194_5e+05_long", 
#"run_Brain_Caudate_basal_ganglia_194_5e+05_short")
dirs = c("run_Brain_Caudate_basal_ganglia_194_5e+05_long", 
"run_Brain_Caudate_basal_ganglia_194_5e+05_short")
fls = list.files(dirs[1], pattern="time")
length(fls)

library(qvalue)
timi = sprintf("%s/%s", dirs[1], fls[1])
resi = read.table(gsub("_time", "_eqtl", timi), as.is=T, header=T)
eqtl5l = data.frame(matrix(NA, nrow=length(fls), ncol=ncol(resi)))
colnames(eqtl5l) = colnames(resi)
eqtl2s = eqtl2l = eqtl5s = eqtl5l
qvals = times = matrix(NA, nrow=length(fls), ncol=4)
for(i in 1:length(fls)){
  #2e5 long
  timi = sprintf("%s/%s", dirs[1], fls[i])
  resi = read.table(gsub("_time", "_eqtl", timi), as.is=T, header=T)
  o = order(resi$final_Pvalue)
  eqtl5l[i,] = resi[o[1],]
  timi = read.table(timi, as.is=T, header=F)
  times[i, 3] = timi[1,1]
  tag = tryCatch({qvals[i, 3] = min(qvalue(na.omit(resi$final_Pvalue))$qvalue);o}, error=function(e){1})

  #2e5 short
  timi = sprintf("%s/%s", dirs[2], fls[i])
  resi = read.table(gsub("_time", "_eqtl", timi), as.is=T, header=T)
  o = order(resi$final_Pvalue)
  eqtl5s[i,] = resi[o[1],]
  timi = read.table(timi, as.is=T, header=F)
  times[i, 4] = timi[1,1]
  tag = tryCatch({qvals[i, 4] = min(qvalue(na.omit(resi$final_Pvalue))$qvalue);o}, error=function(e){1})

  times[i,]
  eqtl2s[i,]
  eqtl2l[i,]
  eqtl5s[i,]
  eqtl5l[i,]
  qvals[i,]
  if(i %% 50 == 0)message(i, " out of ", length(fls))
}
times[1:5,]

fig.dir = "fig"
coli = rgb(.5, .5, .5, .5)
if(!file.exists(fig.dir))dir.create(fig.dir)
pdf(sprintf("%s/times.pdf",fig.dir), height=4, width=4)
par(mfrow=c(1,1))
x = log10(1+times[,4])
y = log10(1+times[,3]) 
plot(x,y, xlab="short", ylab="long", main="5e5 long vs short", bty="n", col=coli, cex=.2)
legend("topleft", legend=c("sec/gene", legend=round(colMeans(times)[c(3,4)])), bty="n")
lm1 = lm(y~x)
abline(a=lm1$coef[1], b=lm1$coef[2])
legend("bottomright", legend=sprintf("%s=%s", c("a", "b"), round(lm1$coef,2)), bty="n")
dev.off()

#eqtl2s[1:5,]
#eqtl2l[1:5,]
#eqtl5s[1:5,]
#eqtl5l[1:5,]

#pdf(sprintf("%s/eqtls.pdf",fig.dir), height=8, width=8)

cutoff = 20
pdf(sprintf("%s/pvals.pdf",fig.dir), height=4, width=4)
par(mfrow=c(1,1))
x = -log10(eqtl5s$final_Pvalue)
y = -log10(eqtl5l$final_Pvalue)
x[x>cutoff] = cutoff
y[y>cutoff] = cutoff
plot(x, y, xlab="short", ylab="long", main="5e5 long vs short", bty="n", col=coli, cex=.2)
abline(a=0, b=1, lty=3, col="red", lwd=2)
dev.off()

cutoff = 5
for(cutq in c(0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6)){
pdf(sprintf("%s/eff_%s.pdf",fig.dir, cutq), height=4, width=12)
par(mfrow=c(1,3))
x = (eqtl5s$TReC_b)
y = (eqtl5l$TReC_b)
x[x>cutoff] = cutoff
y[y>cutoff] = cutoff
x[x< -cutoff] = -cutoff
y[y< -cutoff] = -cutoff
kp = which(qvals[,3]<cutq & qvals[,4]<cutq)
plot(x[kp], y[kp], xlab="short", ylab="long", main="TReC, 5e5 long vs short", bty="n", col=coli, cex=.2)
abline(a=0, b=1, lty=3, col="red", lwd=2)

x = (eqtl5s$ASE_b)
y = (eqtl5l$ASE_b)
x[x>cutoff] = cutoff
y[y>cutoff] = cutoff
x[x< -cutoff] = -cutoff
y[y< -cutoff] = -cutoff
plot(x[kp], y[kp], xlab="short", ylab="long", main="ASE, 5e5 long vs short", bty="n", col=coli, cex=.2)
abline(a=0, b=1, lty=3, col="red", lwd=2)

x = (eqtl5s$Joint_b)
y = (eqtl5l$Joint_b)
x[x>cutoff] = cutoff
y[y>cutoff] = cutoff
x[x< -cutoff] = -cutoff
y[y< -cutoff] = -cutoff
plot(x[kp], y[kp], xlab="short", ylab="long", main="Joint, 5e5 long vs short", bty="n", col=coli, cex=.2)
abline(a=0, b=1, lty=3, col="red", lwd=2)
dev.off()


pdf(sprintf("%s/qvals.pdf",fig.dir), height=4, width=4)
par(mfrow=c(1,1))
x = -log10(qvals[,4])
y = -log10(qvals[,3])
x[x>cutoff] = cutoff
y[y>cutoff] = cutoff
plot(x, y, xlab="short", ylab="long", main="5e5 long vs short", bty="n", col=coli, cex=.2)
abline(a=0, b=1, lty=3, col="red", lwd=2)
dev.off()
}

rownames(qvals) = rownames(pvals) 
write.csv(qvals, "qvals.csv", row.names=T