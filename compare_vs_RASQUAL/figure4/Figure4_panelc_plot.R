workdir = "/pine/scr/z/h/zhabotyn/R01/2019_03_20"
#workdir = "C:/research/R01/reqtl/2019_03_20/sim_ex/snp4/"
setwd(workdir)

res = read.csv("figure4c_summary.csv")
pars = res[,1:5]
pows = res[,6:7]
fule = res[,8:10]
shoe = res[,11:13]


cols = c("black", rep("blue", 3))
kp1 = pars$nsnpi == 1
kp2 = pars$nsnpi == 2
kp3 = pars$nsnpi == 4
x0 = 1:4
x = pars$NBod[kp1]

png("Figure4c.png", height=4, width=4, res=300, units="in")
plot(x0, pows$full[kp1], ylim=c(0,.3), type="b", 
main="(c)BB(OD)=NB(OD)/8", xlab="NB(OD)", bty="n", xaxt="n")
legend("topleft", legend=c("TReCASE", "1 fSNP", "2 fSNPs", "4 fSNPs"), 
bty="n", text.col=cols, col=cols, lty=1:4, pch=1:4,)

axis(1, at=x0, labels=x)

lines(x0, pows$short[kp1], type="b", lty=2, pch=2, col="blue")
lines(x0, pows$short[kp2], type="b", lty=3, pch=3, col="blue")
lines(x0, pows$short[kp3], type="b", lty=4, pch=4, col="blue")
abline(h=0.05, col="red", lty=3)
dev.off()

