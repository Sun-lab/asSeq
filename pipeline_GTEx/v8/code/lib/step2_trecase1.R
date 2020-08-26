args=(commandArgs(TRUE))
#args = c("22", "1", "5e5", "long")
#args = c("12", "496", "467", "5e+05", "long")
chri = as.numeric(args[1])
geni = as.numeric(args[2])
#nsub = as.numeric(args[3])
cis_window = as.numeric(args[3])
model = args[4]
#cis_window = 2e5
chri 
geni
#nsub
cis_window
model

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
asSeq.lib = specs[18]
rcmd = specs[19]
setwd(wrk.dir)

nsub = nsam
c(pref, nsub)

permute = F
set.seed(seedval)

# ------------------------------------------------------------
# read in data
# ------------------------------------------------------------
#root.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8"
#gtex.dir = root.dir

outdir0 = sprintf("run_%s_%s_%s_%s", pref, nsub, cis_window, model)
if(!file.exists(outdir0))dir.create(outdir0)

cnt.dir = sprintf("%s_prepr", pref);cnt.dir;file.exists(cnt.dir)
int.dir = sprintf("%s_%s_%s", pref, nsub, cis_window)


geneInfo = read.table(sprintf("%s/geneInfo_prepr_%s.txt", cnt.dir, model), 
                      header = T, as.is = T)
geneInfo = geneInfo[geneInfo$chr==sprintf("chr%s", chri),]
#geneInfo[1:4,]

library(asSeq, lib.loc=asSeq.lib)


eChr = rep(chri, nrow(geneInfo))
ePos = as.numeric((geneInfo[,3] + geneInfo[,4])/2)


suffi = sprintf("%s_%s", chri, geni)
#genotypes_%s.dat
genfi = read.table(sprintf("%s/genotypes_%s.dat", int.dir, suffi), as.is=T, header=T)
infi = read.table(sprintf("%s/genotypei_%s.dat", int.dir, suffi), as.is=T, header=T)
#genfi[1:4,1:4]
table(unlist(genfi[,1]))
head(infi)

cnts = read.csv(sprintf("%s/counti_%s_%s.csv", int.dir, model, suffi), as.is=T, header=F)
#cnts[1:4,]
Xmatfil = sprintf("%s/Xmat_%s.csv", int.dir, model)
X = read.csv(Xmatfil, as.is=T, header=F)
nr = nrow(genfi)
genfi = matrix(unlist(genfi), nrow=nr)
nsam = nrow(X)
X = matrix(unlist(X), nrow=nsam)
dim(X)
#X[1:4,]

converge = 1e-4
vari = apply(X,2,var)

updvar = which(vari<converge)
for(i in updvar){
  if(length(vari[-updvar]>0)>0){
    correct = sqrt(median(vari[-updvar]))/sqrt(vari[i])
  }else{
    correct = 1/sqrt(vari[i])
  }
  xm = mean(X[,i])    
  X[,i] = xm+(X[,i]-xm)*correct
}
vari
apply(X,2,var)

cnts = matrix(as.numeric(unlist(cnts)), nrow=nsam)
#cnts[1:4,]


mChr = rep(chri, nr)
#geni = 1
#message("geno: ", nrow(geno), " ", ncol(geno), " trecD: ", nrow(trecD), " ", ncol(trecD), " ", nrow(geneInfo))
#geni = 389

#geni = indi
  local.distances = as.numeric((geneInfo[,4] - geneInfo[,3]))/2+cis_window
#  kp = which(SNPInfo[,3]>=(ePos[geni]-local.distances[geni]) & SNPInfo[,3]<(ePos[geni]+local.distances[geni]));length(kp)
  kp = nrow(infi)
  if(kp>0){
    output.tagi     = sprintf("%s/%s",
                               outdir0, geneInfo[geni,1])#rownames(trecD)[geni])
    timj = sprintf("%s_time.txt", output.tagi)
    res.lon = sprintf("%s_eqtl.txt", output.tagi)
      
    time1 = proc.time()
    asSeq:::trecase(Y=cnts[,1,drop=F], Y1=cnts[,2,drop=F], Y2=cnts[,3,drop=F], X=X, Z=t(genfi), output.tag = output.tagi, 
                p.cut=1,  local.distance = local.distances[geni], 
                eChr = eChr[geni], 
                ePos = ePos[geni], 
                mChr = as.numeric(infi[,2]), mPos = as.numeric(infi[,3]),
                maxit = 4000, 
                min.AS.reads = 5, min.AS.sample = 5, min.n.het = 5)

    time2 = proc.time()
    time2 - time1 
  
    write.table(time2[3]-time1[3], timj, row.names=F, col.names=F, quote=F)
    eqtl = read.table(res.lon, header=T, as.is=T)
  
    eqtl[,1] = infi[eqtl[,2],3]; colnames(eqtl)[1] = "Pos"
    eqtl[,2] = infi[eqtl[,2],1]
    write.table(eqtl, res.lon, row.names=F, col.names=T, quote=F, sep="\t")
  }
  message(geni, " out of ", nrow(geneInfo), " #snps: ", length(kp))

gc()

q(save = 'no')


