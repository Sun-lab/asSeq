#sbatch -p general -t 03-00:00:00 -o step6_glm1_1_2.out --mem=4g --wrap="R CMD BATCH '--args 1' step6_glm1.R step6_glm1_1_2.Rout"
args=(commandArgs(TRUE))
# args = c("1", "2")
# args = c("2", "2")
to_run = as.numeric(args[1])
#covi = as.numeric(args[2])
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
out.dir = sprintf("BetaBin_gPC")
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

m = match(rownames(cntL), resL[,1]);m

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
pheno.sub$age = pheno.sub$agen
pheno.sub$tp53 = log10(nrmL[cand[1],]+1)
pheno.sub$ctcf = log10(nrmL[cand[2],]+1)
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

pheno.sub$gPC1 = pheno.sub$gPC1/sd(pheno.sub$gPC1)
pheno.sub$gPC2 = pheno.sub$gPC2/sd(pheno.sub$gPC2)
pheno.sub$PF1 = pheno.sub$PF1/sd(pheno.sub$PF1)
pheno.sub$PF2 = pheno.sub$PF2/sd(pheno.sub$PF2)
pheno.sub$PF3 = pheno.sub$PF3/sd(pheno.sub$PF3)
pheno.sub$PF4 = pheno.sub$PF4/sd(pheno.sub$PF4)
pheno.sub$PF5 = pheno.sub$PF5/sd(pheno.sub$PF5)

pheno.sub$age = pheno.sub$age/sd(pheno.sub$age)
pheno.sub$tp53 = pheno.sub$tp53/sd(pheno.sub$tp53)
pheno.sub$ctcf = pheno.sub$ctcf/sd(pheno.sub$ctcf)
sd(pheno.sub$age)
sd(pheno.sub$tp53)
sd(pheno.sub$ctcf)

#cnd2nm = colnames(pheno.sub)[9+covi];cnd2nm
#pheno.sub$cnd2 = pheno.sub[,cnd2nm]
#pheno.sub[1:2,]

resL1.bb = resS1.bb = resL2.bb = resS2.bb = resL3.bb = resS3.bb = matrix(NA, nrow=nrow(minsL), ncol=23)
resL1.bin = resS1.bin = matrix(NA, nrow=nrow(minsL), ncol=20)
resL1.qp = resS1.qp = resL2.qp = resS2.qp = resL3.qp = resS3.qp = matrix(NA, nrow=nrow(minsL), ncol=20)
min.sam = 15
min.cnt = 10
        #glmL = glm(cbind(yL, nL)~phenoL$sex+phenoL$agen+phenoL$sex*phenoL$agen, family="quasibinomial")#+pheno$DTHHRDY#+Xsho[,6]+Xsho[,7]
        #glmL2 = glm(cbind(yL, nL)~phenoL$SEX+phenoL$agen+phenoL$SEX*phenoL$agen, family="quasibinomial")#+pheno$DTHHRDY#+Xsho[,6]+Xsho[,7]

      i = 0
for(i in 1:nrow(cntL)){#1:100){#
#for(i in 1:500){
#for(i in 1:1000){
#for(i in 2001:nrow(cntL)){
#for(i in 1:100){
#i = i + 1
        #result1 = tryCatch({
        genoL = minsL[i,]
        proceed = TRUE
        if(any(is.na(minsL[i,])))proceed = FALSE
        if(proceed){
        yL = hap1L[i,]
        nL = yL + hap2L[i,]
        yL[genoL==3] = hap2L[i, genoL==3]
  
        kp = which(genoL %in% c(1:3) & nL>=min.cnt);length(kp)
        kpa = which(genoL %in% c(0:4) & nL>=min.cnt);length(kpa)
        #kp = which(genoL %in% c(1,3) & nL>=min.cnt);length(kp)
        if(length(kp)>=min.sam){proceed=TRUE}else{proceed=FALSE}
        }
        proceed
        if(proceed){
          attempt = 1
        table(kp)
        nL = nL[kp]
        yL = yL[kp]
        snpL = genoL[kp]
        #cbind(nL, yL, snpL)
        codL = snpL
        codL[codL==4] = 0
        codL[codL==3] = 1
        phenoL = pheno.sub[kp,]
        
        #phenoL$int = phenoL$sex*phenoL$agen
        phenoL$cod = codL
        #phenoL$cnd = phenoL$agen
        #phenoL$cnd = phenoL$tp53
        phenoL$cnd = phenoL[,conds[to_run]]
        phenoL$int = phenoL$sex*phenoL$cnd
        attempt = tryCatch({
          fit0 = fit1 = fit1a = fit2 = fit = NULL
          fit0 <- vglm(cbind(yL, nL-yL) ~ 1, family=betabinomial, data=phenoL,
              trace = FALSE, subset = nL > 0)

          ini = c(fit0@coefficients, 0, 0); ini[2][ini[2]< -5] = -5
          fit1 <- vglm(cbind(yL, nL-yL) ~ 1 + gPC1 + gPC2, family=betabinomial, data=phenoL,
              trace = FALSE, subset = nL > 0)

          ini = c(fit1@coefficients, 0); ini[2][ini[2]< -5] = -5
          ini[1][ini[1]< -5] = -5
          ini[1][ini[1]> 5] = 5
          fit <- vglm(cbind(yL, nL-yL) ~ 1 + gPC1 + gPC2 + cnd, family=betabinomial, data=phenoL,
              trace = F, subset = nL > 0, coefstart=ini)
#          summary(fit)
          0
        },warning = function(w){
          1
        },error = function(e){
          2
        })
        attempt
        resL1.bb[i,23] = attempt

        if(is.null(fit)){
        tryCatch({
          fit <- vglm(cbind(yL, nL-yL) ~ 1 + gPC1 + gPC2 + cnd, family=betabinomial, data=phenoL,
              trace = F, subset = nL > 0)
        },error = function(e){
          2
        })
        }
        fit.bin1 <- glm(cbind(yL, nL-yL) ~ 1 + gPC1 + gPC2, family=binomial, data=phenoL)
        fit.bin <- glm(cbind(yL, nL-yL) ~ 1 + gPC1 + gPC2 + cnd, family=binomial, data=phenoL)
        fit.qua <- glm(cbind(yL, nL-yL) ~ 1 + gPC1 + gPC2 + cnd, family=quasibinomial, data=phenoL)
        #fit
        #logLik(fit.bin)
        #resL1.bb[i,]
        
        #data.frame(id=rownames(cntL)[i], gni=i, numsam=length(kp), allsam=length(kpa), attempt)
        if(attempt==0){          
          tcnd = 2*(fit@criterion$loglikelihood-fit1@criterion$loglikelihood)
          resL1.bb[i, 1:length(fit@coefficients)] = fit@coefficients
          resL1.bb[i, 10+1:length(summary(fit)@coef3[,4])] = summary(fit)@coef3[,4]
          resL1.bb[i, 21] = pchisq(tcnd, df=1, lower.tail=F)
        }
        if(!is.null(fit)){
          tryCatch({        
            resL1.bb[i, 22] = logLik(fit)
            resL1.bb[i, 10+1:length(summary(fit)@coef3[,4])] = summary(fit)@coef3[,4]
            resL1.bb[i, 1:length(fit@coefficients)] = fit@coefficients
            1
          },error = function(e){
            2
          })        
        }
        resL1.bin[i, 1:length(fit.bin$coefficients)] = fit.bin$coefficients
        resL1.bin[i, 9+1:nrow(summary(fit.bin)$coefficients)] = summary(fit.bin)$coefficients[,4]
        resL1.bin[i, 20] = logLik(fit.bin)

          tcnd = 2*(logLik(fit.bin)-logLik(fit.bin1))
          resL1.bin[i, 19] = pchisq(tcnd, df=1, lower.tail=F)


        resL1.qp[i, 1:length(fit.qua$coefficients)] = fit.qua$coefficients
        resL1.qp[i, 9+1:nrow(summary(fit.qua)$coefficients)] = summary(fit.qua)$coefficients[,4]
        resL1.qp[i, 20] = summary(fit.qua)$dispersion
        }

        genoS = minsS[i,]
        proceed = TRUE
        if(any(is.na(minsS[i,])))proceed = FALSE
        if(proceed){
        yS = hap1S[i,]
        nS = yS + hap2S[i,]
        yS[genoS==3] = hap2S[i, genoS==3]
  
        kp = which(genoS %in% c(1:3) & nS>=min.cnt);length(kp)
        kpa = which(genoS %in% c(0:4) & nS>=min.cnt);length(kpa)
        #kp = which(genoS %in% c(1,3) & nS>=min.cnt);length(kp)
        if(length(kp)>=min.sam){proceed=TRUE}else{proceed=FALSE}
        }
        proceed
        if(proceed){
          attempt = 1
        table(kp)
        nS = nS[kp]
        yS = yS[kp]
        snpS = genoS[kp]
        #cbind(nS, yS, snpS)
        codS = snpS
        codS[codS==4] = 0
        codS[codS==3] = 1
        phenoS = pheno.sub[kp,]
        #phenoS$int = phenoS$sex*phenoS$agen
        phenoS$cod = codS
        #phenoS$cnd = phenoS$agen
        #phenoS$cnd = phenoS$tp53
        phenoS$cnd = phenoS[,conds[to_run]]
        phenoS$int = phenoS$sex*phenoS$cnd
        attempt = tryCatch({
          fit0 = fit1 = fit1a = fit2 = fit = NULL

          fit0 <- vglm(cbind(yS, nS-yS) ~ 1, family=betabinomial, data=phenoS,
              trace = FALSE, subset = nS > 0)

          ini = c(fit0@coefficients, 0, 0); ini[2][ini[2]< -5] = -5
          fit1 <- vglm(cbind(yS, nS-yS) ~ 1 + gPC1 + gPC2, family=betabinomial, data=phenoS,
              trace = FALSE, subset = nS > 0)

          ini = c(fit1@coefficients, 0); ini[2][ini[2]< -5] = -5
          ini[1][ini[1]< -5] = -5
          ini[1][ini[1]> 5] = 5
          fit <- vglm(cbind(yS, nS-yS) ~ 1 + gPC1 + gPC2 + cnd, family=betabinomial, data=phenoS,
              trace = F, subset = nS > 0, coefstart=ini)
          0
        },warning = function(w){
          1
        },error = function(e){
          2
        })
        attempt
        resS1.bb[i,23] = attempt

        if(is.null(fit)){
        tryCatch({
          fit <- vglm(cbind(yS, nS-yS) ~ 1 + gPC1 + gPC2 + cnd, family=betabinomial, data=phenoS,
              trace = F, subset = nS > 0)
        },error = function(e){
          2
        })
        }
        fit.bin1 <- glm(cbind(yS, nS-yS) ~ 1 + gPC1 + gPC2, family=binomial, data=phenoS)
        fit.bin <- glm(cbind(yS, nS-yS) ~ 1 + gPC1 + gPC2 + cnd, family=binomial, data=phenoS)
        fit.qua <- glm(cbind(yS, nS-yS) ~ 1 + gPC1 + gPC2 + cnd, family=quasibinomial, data=phenoS)
        #fit
        #logLik(fit.bin)
        #resS1.bb[i,]
        
        #data.frame(id=rownames(cntS)[i], gni=i, numsam=length(kp), allsam=length(kpa), attempt)
        if(attempt==0){          
          tcnd = 2*(fit@criterion$loglikelihood-fit1@criterion$loglikelihood)
          resS1.bb[i, 1:length(fit@coefficients)] = fit@coefficients
          resS1.bb[i, 10+1:length(summary(fit)@coef3[,4])] = summary(fit)@coef3[,4]
          resS1.bb[i, 21] = pchisq(tcnd, df=1, lower.tail=F)
        }
        if(!is.null(fit)){
          tryCatch({        
            resS1.bb[i, 22] = logLik(fit)
            resS1.bb[i, 1:length(fit@coefficients)] = fit@coefficients
            resS1.bb[i, 10+1:length(summary(fit)@coef3[,4])] = summary(fit)@coef3[,4]
            1
          },error = function(e){
            2
          })        
        }
        resS1.bin[i, 1:length(fit.bin$coefficients)] = fit.bin$coefficients
        resS1.bin[i, 9+1:nrow(summary(fit.bin)$coefficients)] = summary(fit.bin)$coefficients[,4]
        resS1.bin[i, 20] = logLik(fit.bin)

          tcnd = 2*(logLik(fit.bin)-logLik(fit.bin1))
          resS1.bin[i, 19] = pchisq(tcnd, df=1, lower.tail=F)


        resS1.qp[i, 1:length(fit.qua$coefficients)] = fit.qua$coefficients
        resS1.qp[i, 9+1:nrow(summary(fit.qua)$coefficients)] = summary(fit.qua)$coefficients[,4]
        resS1.qp[i, 20] = summary(fit.qua)$dispersion
        }
        if(i%%100==0)message(i, " out of ", nrow(cntS))      
  }

rownames(resL1.qp) = rownames(resS1.qp) =
rownames(resL1.bin) = rownames(resS1.bin) =
rownames(resL1.bb) = rownames(resL2.bb) = rownames(resL3.bb) = 
rownames(resS1.bb) = rownames(resS2.bb) = rownames(resS3.bb) = rownames(cntL)


mod = sprintf("%s_%s", pref, conds[to_run]);mod
write.csv(resL1.bb, sprintf("%s/%s_long_vglm_gPC.csv", out.dir, mod), quote=F)
write.csv(resS1.bb, sprintf("%s/%s_short_vglm_gPC.csv", out.dir, mod), quote=F)

write.csv(resL1.qp, sprintf("%s/%s_long_glmqp_gPC.csv", out.dir, mod), quote=F)
write.csv(resS1.qp, sprintf("%s/%s_short_glmqp_gPC.csv", out.dir, mod), quote=F)

write.csv(resL1.bin, sprintf("%s/%s_long_glmbin_gPC.csv", out.dir, mod), quote=F)
write.csv(resS1.bin, sprintf("%s/%s_short_glmbin_gPC.csv", out.dir, mod), quote=F)

gc()

q("no")
