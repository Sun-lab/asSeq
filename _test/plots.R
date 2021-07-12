plotASE <- 
function(y1, y2, xx1, xlab, ylab, sub="", min.nTotal=10)
{
  yT = y1 + y2
  y2[xx1==3] = y1[xx1==3]
  
  w2use = which(yT>=min.nTotal)
  y2 = as.numeric(y2[w2use])
  yT = as.numeric(yT[w2use])
  x1 = as.numeric(xx1[w2use])
  
  xh = rep("homozyg", length(x1))
  xh[which(abs(x1-2)<1.5)] = "heterozyg"
  zeta = as.numeric(xh == "heterozyg")
  
  a1 = aseR(y2, yT, zeta, maxIt=50, trace=0)
  p1 = 1 - pchisq(2.0*a1$logLikH1 - 2.0*a1$logLikH0, 1)
  mm1 = sprintf("pASE=%.1e", p1)
  
  boxplot(y2/yT ~ xh, xlab=xlab, ylab=ylab, ylim=c(0,1), 
          outline=FALSE, main=mm1, bty="n")
  stripchart(y2/yT ~ xh, method = "jitter", jitter = 0.2, offset = 1/3,
             vertical = TRUE, add = TRUE, col = "darkgreen", pch=20)
  abline(h=0.5, col="skyblue", lty=2)
  mtext(sub, side=1, line=2, cex=0.6)
}

plotTReC <- 
function(y0, xx1, X, xlab, ylab)
{  
  y0 = as.numeric(y0)

  z1 = as.numeric(xx1)
  z1[which(xx1>=3)] = z1[which(xx1>=3)] - 2
    
  X0 = rep(1, length(y0))
  l0 = trecR(y0, X0, z1, fam="negbin")
  l1 = trecR(y0, X,  z1, fam="negbin", yfit=TRUE)
  
  p0 = 1 - pchisq(l0$lrt, 1)
  p1 = 1 - pchisq(l1$lrt, 1)
  
  boxplot(y0 ~ z1, main=sprintf("Before Adjust. (p=%.1e)", p0), 
          xlab=xlab, ylab=ylab, outline=FALSE, bty="n")
  
  stripchart(y0 ~ z1, method = "jitter", jitter = 0.2, offset = 1/3,
             vertical = TRUE, add = TRUE, col = "darkgreen", pch=20)
  
  y1 = l1$fitted
  
  boxplot(y1 ~ z1, main=sprintf("After Adjust. (p=%.1e)", p1),
          xlab=xlab, ylab=ylab, outline=FALSE)
  
  stripchart(y1 ~ z1, method = "jitter", jitter = 0.2, offset = 1/3,
             vertical = TRUE, add = TRUE, col = "darkgreen", pch=20)
  
}
