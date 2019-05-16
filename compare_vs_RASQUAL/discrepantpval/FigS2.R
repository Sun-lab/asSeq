png("FigS2_upd.png", height=4, width=12, res=300, units="in")
cexes = 1.8
par(mfrow=c(1,3))
par(mar=c(5,5,3,1))

plot(ind, frR[1,], xlab="|log10(T p-val)-log10(R p-val)|", ylab="%p(R)<p(T)", 
main="", bty="n", cex.main=cexes, cex.lab=cexes, cex.axis=cexes, cex=log10(frR[2,]))
legend("bottomright", legend=c("#genes", 10^(4:2)), pt.cex=c(1, 4:2), 
       cex=cexes, bty="n", pch=c(NA, 1,1,1))
legend("topleft", legend=c(2500, 500, 100), pt.cex=log10(c(2500, 500, 100)), 
       cex=cexes, bty="n", pch=c(1,1,1))
axis(1, at=c(12), 12, cex.axis=cexes)


x = 1:6
plot(x, ag0[,2], ylim=c(0,1), type="b", xlab="#fSNP", ylab="%p(R)<p(T)", 
      main="", bty="n", xaxt="n", cex=log10(cg0[,2]), 
      xlim=c(1,6.1), cex.main=cexes, cex.lab=cexes, cex.axis=cexes)
points(x, ag5[,2], col="blue", type="b", pch=2, cex=log10(cg5[,2]))
points(x, ag10[,2], col="goldenrod", type="b", pch=3, cex=log10(cg10[,2]))
points(x, ag15[,2], col="red", type="b", pch=4, cex=log10(cg15[,2]))
axis(1, at=x, ag0[,1], cex.lab=cexes, cex.axis=cexes)
legend("topleft", legend=c("all",">5", ">10", ">15"), 
text.col=c("black", "blue", "goldenrod", "red"), pch=c(1:4), bty="n", cex=cexes)


plot(ind, cors, cex=log10(numg), pch=1, bty="n", xlab="|log10(T p-val)-log10(R p-val)|", 
ylab="OD correlation", main="", cex.main=cexes, cex.lab=cexes, cex.axis=cexes, ylim=c(0,1))
points(ind, cors2, cex=log10(numg), pch=2)
axis(1, at=c(12), 12, cex.axis=cexes)
legend("bottomleft", c("OD(NB)", "OD(BB)"), bty="n", pch=c(1,2), cex=cexes)
legend("bottomright", legend=c("#genes", 10^(4:2)), 
       pt.cex=c(1, 4:2), cex=cexes, bty="n", pch=c(NA, 1,1,1))
dev.off()

