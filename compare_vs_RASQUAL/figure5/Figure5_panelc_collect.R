#if you want to run multiple b0's:
#about 8.5 seconds per iteration for all profiles
#expect 71 minutes per 500 iterations
#if you run only b0=0 - under 15 minutes

niter = 5e2
queue = "general"
mem = 8000
ss = 2
#b0s = c(0, .125,.25,.5, 1)
b0s = c(0, .5, 1)
#div = c(.5, 1, 2, 4, 8, Inf)
div = 1#a subset for figure (a)
ods = c(0.01, 0.10, 0.5, 2)

days = 7
mem = 8
b1 = 0
nsnp = c(4)
nsim = length(ods)*length(div)*length(b0s)*length(ss)*length(nsnp)
pars = matrix(NA, nrow=nsim, ncol=5)
fule = matrix(NA, nrow=nsim, ncol=3)
shoe = matrix(NA, nrow=nsim, ncol=3)
pows = matrix(NA, nrow=nsim, ncol=2)
colnames(pows) = c("full", "short")
colnames(pars) = c("sample.s", "b0", "divisor", "NBod", "nsnpi")
colnames(fule) = c("fullNBod", "fullBBod", "fullb0")
colnames(shoe) = c("shrtNBod", "shrtBBod", "shrtb0")
i = 0

alpha = 0.05
powerpos = 2+1+4+1+1#2 ods, 1 b0, 4 betas, 1 ll, p-val 
st = proc.time()
#nsnpi = 1; ssi = 2; b0i = 0; dvi = 1; odi = 0.5
for(nsnpi in nsnp){
  for(ssi in ss){
  for(b0i in b0s){
      for(dvi in div){
        for(odi in ods){
          i = i + 1
          pars[i,] = c(ssi, b0i, dvi, odi, nsnpi)
          suff = sprintf("add%s_phi%s_div%s_iter%s_ss%s_nsnp%s", 
                         b0i, odi, dvi, niter, ssi, nsnpi)
          fm2od = read.csv(sprintf("resnew/full_%s.csv", suff), header=F)
          fm1od = read.csv(sprintf("resnew/short_%s.csv", suff),header=F)

          pows[i,] = c(mean(fm2od[,powerpos]<alpha), mean(fm1od[,powerpos]<alpha))
          fule[i,] = apply(fm2od[,1:3], 2, median)
          shoe[i,] = apply(fm1od[,1:3], 2, median)
          fule[i,1] = -fule[i,1]
          shoe[i,1] = -shoe[i,1]
          fule[i,1:2] = exp(fule[i,1:2])          
          shoe[i,1:2] = exp(shoe[i,1:2])
          message(ssi, " ", b0i, " ", dvi, " ", odi, " ", nsnpi)      
          }
        }
      }
  }
}
en = proc.time()
en[3]-st[3]

write.csv(cbind(pars, pows, fule, shoe), "figure5c_summary.csv", row.names=F)
