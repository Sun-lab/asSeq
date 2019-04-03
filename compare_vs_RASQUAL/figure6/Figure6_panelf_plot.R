
#workdir = "/pine/scr/z/h/zhabotyn/R01/2019_03_20"
#workdir = "C:/Users/Vasyl/Documents/GitHub/asSeq/compare_vs_RASQUAL/figure6";setwd(workdir)

res = read.csv("figure6f_summary.csv")
pars = res[,1:5]
pows = res[,6:7]
fule = res[,8:10]
shoe = res[,11:13]


ods = c(0.01, 0.1, 0.5, 2)
cols = c("black", rep("blue", 3))
i = 1
kp = pars$NBod == ods[i]
x = pars$b0[kp]

png("Figure6f.png", height=4, width=4, res=300, units="in")
plot(x, fule$fullb0[kp], ylim=c(-0.05,1.05), xlim=c(-0.05,1.05), type="b", 
main="(f)eQTL, 4SNPs, BB(OD)=0", xlab="true eQTL", ylab="median(est eQTL)", bty="n")
lines(x, shoe$shrtb0[kp], type="b", pch=2, col="blue")
legend("bottomright", legend=c("TReCASE", "TReCASE-RL"), 
bty="n", text.col=cols, col=cols, lty=1:2, pch=1:2)
legend("topleft", legend=c("OD", ods), bty="n", pch=1, , pt.cex=1+c(0,0,(2:4)/4),col=c("white", rep("black", 4)))
for(i in 2:4){
  kp = pars$NBod == ods[i]
  lines(x, fule$fullb0[kp], type="b", lty=1, pch=1, cex=1+i/4, col="black")
  lines(x, shoe$shrtb0[kp], type="b", lty=2, pch=2, cex=1+i/4, col="blue")
}
abline(a=0, b=1, col="red", lwd=3)
dev.off()


