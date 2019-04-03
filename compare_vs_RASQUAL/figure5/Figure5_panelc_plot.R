#workdir = "/pine/scr/z/h/zhabotyn/R01/2019_03_20"
#workdir = "C:/research/R01/reqtl/2019_03_20/sim_ex/snp4/"
#setwd(workdir)

res = read.csv("figure5c_summary.csv")
pars = res[,1:5]
pows = res[,6:7]
fule = res[,8:10]
shoe = res[,11:13]


cols = c("black", rep("blue", 3))
kp1 = pars$b0 == 0
kp2 = pars$b0 == .5
kp3 = pars$b0 == 1
x = pars$NBod[kp1]

png("Figure5c.png", height=4, width=4, res=300, units="in")
plot(x, fule$fullNB[kp1], type="b", 
main="(c)BB(OD)=NB(OD)", xlab="NB(OD)", bty="n")
legend("topleft", legend=c("TReCASE", "1 fSNP", "2 fSNPs", "4 fSNPs"), 
bty="n", text.col=cols, col=cols, lty=1:4, pch=1:4,)
lines(x, fule$fullNB[kp2], type="b", lty=1, pch=1, col="black", cex=.5+pars$b0[kp2])
lines(x, fule$fullNB[kp3], type="b", lty=1, pch=1, col="black",cex=.5+pars$b0[kp3])

#axis(1, at=x0, labels=x)

lines(x, shoe$shrtNB[kp1], type="b", lty=2, pch=2, col="blue", cex=.5+pars$b0[kp1])
lines(x, shoe$shrtNB[kp2], type="b", lty=2, pch=2, col="blue", cex=.5+pars$b0[kp2])
lines(x, shoe$shrtNB[kp3], type="b", lty=2, pch=2, col="blue",cex=.5+pars$b0[kp3])
abline(a=0, b=1, col="red", lty=3)
dev.off()

