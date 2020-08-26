#sbatch -p general -t 03-00:00:00 -o step6_glm1_1_2.out --mem=4g --wrap="R CMD BATCH '--args 1' step6_glm1.R step6_glm1_1_2.Rout"
args=(commandArgs(TRUE))
# args = c("1", "2")
# args = c("2", "2")
covi = as.numeric(args[1])
conds = c("age", "tp53", "ctcf")

library(VGAM)
#library("HRQoL")

specf = "specifications.txt"
getwd()

specs = unlist(read.table(specf, as.is=T))
specs
pref = specs[1]
nsam = specs[2]
queue = specs[3]
days = specs[4]
bmem = as.numeric(specs[5])
seedval = specs[13]
wrk.dir = specs[14]
lib.dir = specs[15]
bas.dir = specs[16]
setwd(wrk.dir)

source(sprintf("%s/helpers.R", lib.dir))
#v8.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8"
#tis.dir = sprintf("%s/2020_01_20/%s", v8.dir, pref)
v8.dir = bas.dir
tis.dir = wrk.dir
out.dir = sprintf("%s/BetaBin_marginal", wrk.dir)
if(!file.exists(out.dir))dir.create(out.dir)
#setwd(wrk.dir)

ann.dir = sprintf("%s/Annotations", v8.dir)
cov.dir = sprintf("%s/GTEx_Analysis_v8_eQTL_covariates", ann.dir)
cov.fil = sprintf("%s/%s.v8.covariates.txt",cov.dir, pref)
covars = read.table(cov.fil, header=T, as.is=T)
rownames(covars) = covars$ID; covars = covars[,-1]
colnames(covars) = gsub("\\.", "-", colnames(covars))
covars[1:2,1:5]

#cov.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8/Annotations/GTEx_Analysis_v8_eQTL_covariates"
#covs = read.table(sprintf("%s/Brain_Caudate_basal_ganglia.v8.covariates.txt", cov.dir), as.is=T)
phen.fil = sprintf("%s/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", ann.dir)
pheno = read.table(phen.fil, as.is=T, sep="\t", header=T)
pheno[1:4,]

pheno$agen = 25
pheno$agen[pheno$AGE=="20-29"] = 25
pheno$agen[pheno$AGE=="30-39"] = 35
pheno$agen[pheno$AGE=="40-49"] = 45
pheno$agen[pheno$AGE=="50-59"] = 55
pheno$agen[pheno$AGE=="60-69"] = 65
pheno$agen[pheno$AGE=="70-79"] = 75
pheno$sex = -1
pheno$sex[pheno$SEX==2] = 1

#reduce redundancy - in previous step already produced needed conversion of counts
#these lines are still needed though - if we include for fitting only permp<cutoff
#lonf = sprintf("%s/short_TReCASE_%s_%s_5e+05_long.csv", tis.dir, pref, nsam)
#shof = sprintf("%s/short_TReCASE_%s_%s_5e+05_short.csv", tis.dir, pref, nsam)
lonf = sprintf("%s/%s_%s_5e+05_long_TReCASE_trimmed.csv", tis.dir, pref, nsam)
shof = sprintf("%s/%s_%s_5e+05_long_TReCASE_trimmed.csv", tis.dir, pref, nsam)
resL = read.csv(lonf, as.is=T)
resS = read.csv(shof, as.is=T)
resL[1:4,]
resS[1:4,]

lonf = sprintf("%s/counts_with_min_snp_%s_%s_5e+05_long.csv", tis.dir, pref, nsam)
shof = sprintf("%s/counts_with_min_snp_%s_%s_5e+05_short.csv", tis.dir, pref, nsam)
cntL = read.csv(lonf, as.is=T)
cntS = read.csv(shof, as.is=T)
colnames(cntL) = gsub("\\.", "-", colnames(cntL))
colnames(cntS) = gsub("\\.", "-", colnames(cntS))
rownames(cntL) = cntL[,1]; cntL = cntL[,-1]; cntL = as.matrix(cntL)
rownames(cntS) = cntS[,1]; cntS = cntS[,-1]; cntS = as.matrix(cntS)
cntL[1:4,1:4]
cntS[1:4,1:4]

cand = c(grep("ENSG00000141510", rownames(cntL)), grep("ENSG00000102974", rownames(cntL)))

minsL = apply(cntL, 1:2, splitting)
minsS = apply(cntS, 1:2, splitting)

totL = apply(cntL, 1:2, splitting, block=1)
crv = colSums(totL); crv = median(crv)/crv
nrmL = totL%*%diag(crv);summary(colSums(nrmL))

hap1L = apply(cntL, 1:2, splitting, block=2)
hap2L = apply(cntL, 1:2, splitting, block=3)
hap1S = apply(cntS, 1:2, splitting, block=2)
hap2S = apply(cntS, 1:2, splitting, block=3)

sam.ord = colnames(cntL)
m = match(sam.ord, pheno$SUBJID)
table(pheno$SUBJID[m] == sam.ord)
pheno.sub = pheno[m,]
#pheno.sub$age = pheno.sub$agen
#pheno.sub$tp53 = log10(nrmL[cand[1],]+1)
#pheno.sub$ctcf = log10(nrmL[cand[2],]+1)
pheno.sub$age = NA
pheno.sub$tp53 = NA
pheno.sub$ctcf = NA
cutoff = 0.01

m = match(sam.ord, colnames(covars))
covars = covars[,m]
table(sam.ord==colnames(covars))
pheno.sub$gPC1 = unlist(covars["PC1",])
pheno.sub$gPC2 = unlist(covars["PC2",])
pheno.sub$PF1 = unlist(covars["InferredCov1",])
pheno.sub$PF2 = unlist(covars["InferredCov2",])
pheno.sub$PF3 = unlist(covars["InferredCov3",])
pheno.sub$PF4 = unlist(covars["InferredCov4",])
pheno.sub$PF5 = unlist(covars["InferredCov5",])
pheno.sub$PF6 = unlist(covars["InferredCov6",])
pheno.sub$PF7 = unlist(covars["InferredCov7",])
pheno.sub$PF8 = unlist(covars["InferredCov8",])
pheno.sub$PF9 = unlist(covars["InferredCov9",])
pheno.sub$PF10 = unlist(covars["InferredCov10",])

pheno.sub$gPC1 = pheno.sub$gPC1/sd(pheno.sub$gPC1)
pheno.sub$gPC2 = pheno.sub$gPC2/sd(pheno.sub$gPC2)
pheno.sub$PF1 = pheno.sub$PF1/sd(pheno.sub$PF1)
pheno.sub$PF2 = pheno.sub$PF2/sd(pheno.sub$PF2)
pheno.sub$PF3 = pheno.sub$PF3/sd(pheno.sub$PF3)
pheno.sub$PF4 = pheno.sub$PF4/sd(pheno.sub$PF4)
pheno.sub$PF5 = pheno.sub$PF5/sd(pheno.sub$PF5)
pheno.sub$PF6 = pheno.sub$PF6/sd(pheno.sub$PF6)
pheno.sub$PF7 = pheno.sub$PF7/sd(pheno.sub$PF7)
pheno.sub$PF8 = pheno.sub$PF8/sd(pheno.sub$PF8)
pheno.sub$PF9 = pheno.sub$PF9/sd(pheno.sub$PF9)
pheno.sub$PF10 = pheno.sub$PF10/sd(pheno.sub$PF1)

cnd2nm = colnames(pheno.sub)[9+covi];cnd2nm
pheno.sub$cnd2 = pheno.sub[,cnd2nm]
pheno.sub[1:2,]

resL1.bb = resS1.bb = matrix(NA, nrow=nrow(minsL), ncol=9)
resL1.bin = resS1.bin = matrix(NA, nrow=nrow(minsL), ncol=6)
min.sam = 10
min.cnt = 10

i = 0
for(i in 1:nrow(cntL)){#1:100){#
#for(i in 1:100){
#i = i + 1
      genoL = minsL[i,]
      proceed = TRUE
      if(any(is.na(minsL[i,])))proceed = FALSE
      if(proceed){
        yL = hap1L[i,]
        nL = yL + hap2L[i,]
        yL[genoL==3] = hap2L[i, genoL==3]
  
        kp = which(genoL %in% c(1:3) & nL>=min.cnt);length(kp)
        kpa = which(genoL %in% c(0:4) & nL>=min.cnt);length(kpa)
        if(length(kp)>=min.sam){proceed=TRUE}else{proceed=FALSE}
      }
      proceed
      if(proceed){
          attempt = 1
        table(kp)
        nL = nL[kp]
        yL = yL[kp]
        snpL = genoL[kp]
        codL = snpL
        codL[codL==4] = 0
        codL[codL==3] = 1
        phenoL = pheno.sub[kp,]
        
        phenoL$cod = codL
#        phenoL$cnd = phenoL[,conds[to_run]]
#        phenoL$int = phenoL$sex*phenoL$cnd
        attempt = tryCatch({
          fit0 = fit1 = fit1a = fit2 = fit = NULL
          fit0 <- vglm(cbind(yL, nL-yL) ~ 1, family=betabinomial, data=phenoL,
              trace = FALSE, subset = nL > 0)

          ini = c(fit0@coefficients, 0); ini[2][ini[2]< -5] = -5
          ini[1][ini[1]< -5] = -5
          ini[1][ini[1]> 5] = 5
          fit1 <- vglm(cbind(yL, nL-yL) ~ 1+cnd2, family=betabinomial, data=phenoL,
              trace = FALSE, subset = nL > 0)
          0
        },warning = function(w){
          1
        },error = function(e){
          2
        })
        attempt
        resL1.bb[i,9] = attempt

        if(is.null(fit1)){
        tryCatch({
          fit1 <- vglm(cbind(yL, nL-yL) ~ 1+cnd2, family=betabinomial, data=phenoL,
              trace = FALSE, subset = nL > 0)
        },error = function(e){
          2
        })
        }
        fit.bin0 <- glm(cbind(yL, nL-yL) ~ 1, family=binomial, data=phenoL)
        fit.bin1 <- glm(cbind(yL, nL-yL) ~ 1+ cnd2, family=binomial, data=phenoL)
        if(attempt==0){          
          teqt = 2*(fit1@criterion$loglikelihood-fit0@criterion$loglikelihood)# other cov
          resL1.bb[i, 1:length(fit1@coefficients)] = fit1@coefficients
          resL1.bb[i, 3+1:nrow(summary(fit1)@coef3)] = summary(fit1)@coef3[,4]
          resL1.bb[i, 7] = pchisq(teqt, df=1, lower.tail=F)
        }
        if(!is.null(fit1)){
          tryCatch({        
            resL1.bb[i, 3+1:nrow(summary(fit1)@coef3)] = summary(fit1)@coef3[,4]
            resL1.bb[i, 1:length(fit1@coefficients)] = fit1@coefficients
            resL1.bb[i, 8] = logLik(fit1)
            1
          },error = function(e){
            2
          })        
        }
        resL1.bin[i, 1:length(fit.bin1$coefficients)] = fit.bin1$coefficients
        resL1.bin[i, 2+1:nrow(summary(fit.bin1)$coefficients)] = summary(fit.bin1)$coefficients[,4]
        resL1.bin[i, 6] = logLik(fit.bin1)

        teqt = 2*(logLik(fit.bin1)-logLik(fit.bin0))
        resL1.bin[i, 5] = pchisq(teqt, df=1, lower.tail=F)
      }

      #short
      genoS = minsS[i,]
      proceed = TRUE
      if(any(is.na(minsS[i,])))proceed = FALSE
      if(proceed){
        yS = hap1S[i,]
        nS = yS + hap2S[i,]
        yS[genoS==3] = hap2S[i, genoS==3]
  
        kp = which(genoS %in% c(1:3) & nS>=min.cnt);length(kp)
        kpa = which(genoS %in% c(0:4) & nS>=min.cnt);length(kpa)
        if(length(kp)>=min.sam){proceed=TRUE}else{proceed=FALSE}
      }
      proceed
      if(proceed){
          attempt = 1
        table(kp)
        nS = nS[kp]
        yS = yS[kp]
        snpS = genoS[kp]
        codS = snpS
        codS[codS==4] = 0
        codS[codS==3] = 1
        phenoS = pheno.sub[kp,]
        
        phenoS$cod = codS
#        phenoS$cnd = phenoS[,conds[to_run]]
#        phenoS$int = phenoS$sex*phenoS$cnd
        attempt = tryCatch({
          fit0 = fit1 = fit1a = fit2 = fit = NULL
          fit0 <- vglm(cbind(yS, nS-yS) ~ 1, family=betabinomial, data=phenoS,
              trace = FALSE, subset = nS > 0)

          ini = c(fit0@coefficients, 0); ini[2][ini[2]< -5] = -5
          ini[1][ini[1]< -5] = -5
          ini[1][ini[1]> 5] = 5
          fit1 <- vglm(cbind(yS, nS-yS) ~ 1+cnd2, family=betabinomial, data=phenoS,
              trace = FALSE, subset = nS > 0)
          0
        },warning = function(w){
          1
        },error = function(e){
          2
        })
        attempt
        resS1.bb[i,9] = attempt

        if(is.null(fit1)){
        tryCatch({
          fit1 <- vglm(cbind(yS, nS-yS) ~ 1+cnd2, family=betabinomial, data=phenoS,
              trace = FALSE, subset = nS > 0)
        },error = function(e){
          2
        })
        }
        fit.bin0 <- glm(cbind(yS, nS-yS) ~ 1, family=binomial, data=phenoS)
        fit.bin1 <- glm(cbind(yS, nS-yS) ~ 1+ cnd2, family=binomial, data=phenoS)

        if(attempt==0){          
          teqt = 2*(fit1@criterion$loglikelihood-fit0@criterion$loglikelihood)# other cov
          resS1.bb[i, 1:length(fit1@coefficients)] = fit1@coefficients
          resS1.bb[i, 3+1:nrow(summary(fit1)@coef3)] = summary(fit1)@coef3[,4]
          resS1.bb[i, 7] = pchisq(teqt, df=1, lower.tail=F)
        }
        if(!is.null(fit1)){
          tryCatch({        
            resS1.bb[i, 3+1:nrow(summary(fit1)@coef3)] = summary(fit1)@coef3[,4]
            resS1.bb[i, 1:length(fit1@coefficients)] = fit1@coefficients
            resS1.bb[i, 8] = logLik(fit1)
            1
          },error = function(e){
            2
          })        
        }
        resS1.bin[i, 1:length(fit.bin1$coefficients)] = fit.bin1$coefficients
        resS1.bin[i, 2+1:nrow(summary(fit.bin1)$coefficients)] = summary(fit.bin1)$coefficients[,4]
        resS1.bin[i, 6] = logLik(fit.bin1)

        teqt = 2*(logLik(fit.bin1)-logLik(fit.bin0))
        resS1.bin[i, 5] = pchisq(teqt, df=1, lower.tail=F)
      }
      if(i%%100==0)message(i, " out of ", nrow(cntS))      
}

rownames(resL1.bin) = rownames(resS1.bin) =
rownames(resL1.bb) = rownames(resS1.bb) = rownames(cntL)

colnames(resL1.bb) = colnames(resS1.bb) = 
c("b.i1", "b.i2", "b.cnd2", 
"p.i1", "p.i2", "p.cnd2", 
"pll.cnd2", "logLik", "status")

colnames(resL1.bin) = colnames(resS1.bin) = 
c("b.i1", "b.cnd2", 
"p.i1", "p.cnd2", "pll.cnd2", "logLik")


table(resL1.bb[,"status"])
table(resS1.bb[,"status"])

mod = sprintf("%s_%s", pref, cnd2nm);mod
write.csv(resL1.bb, sprintf("%s/%s_long_vglm_mar.csv", out.dir, mod), quote=F)
write.csv(resS1.bb, sprintf("%s/%s_short_vglm_mar.csv", out.dir, mod), quote=F)

write.csv(resL1.bin, sprintf("%s/%s_long_glmbin_mar.csv", out.dir, mod), quote=F)
write.csv(resS1.bin, sprintf("%s/%s_short_glmbin_mar.csv", out.dir, mod), quote=F)

gc()

q("no")
