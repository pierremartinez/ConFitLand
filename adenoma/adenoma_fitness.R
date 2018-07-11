model="combined_prod"

load(paste("../models/bestCombs_",model,".RData", sep=""))
load(paste("../models/bestFitNL_",model,".RData", sep=""))
load(paste("../models/bestFitL_",model,".RData", sep=""))
## Most clonally mutated genes
load("../tcga/topL.RData")
## list of ratio of mean observed n partners compared to mean expected n partners
load("../tcga/socRatL.RData")
## list of power values for the observed ratios
load("../tcga/socPowL.RData")
## list of selection values 
load("../tcga/selectMutL.RData")

## Data on select CN events
load("../tcga/cnSelectL.RData")

## util functions
getEpi <- function(ints, g, subpanel) {
  partners <- setdiff(subpanel,g)
  if (length(partners) > 0) {
    if (length(grep("_mean$", model)) > 0) {
      return(mean(unlist(ints[g,setdiff(subpanel,g)])))
    } else {
      return(prod(ints[g,setdiff(subpanel,g)]))
    }
  } else {
    return(1)
  }
}

## Combined model function
getFitStr <- function(panelStr, sel, socRat, ints, sf, ef) {
  if(panelStr=="")
    return(0)
  subpanel=unlist(strsplit(panelStr, ":"))
  return(sum(unlist(lapply(subpanel, function(g){sel[g] * socRat[g] * sf + getEpi(ints, g, subpanel) * ef}))))
}

## Update if using separate model
## getFitStr <- function(panelStr, sel, socRat, ints, of, sf, ef) {
##   subpanel=unlist(strsplit(panelStr, ":"))
##   return(sum(unlist(lapply(subpanel, function(g){(sel[g] * of) + (socRat[g] * sf) + (getEpi(ints, g, subpanel) * ef)}))))
## }

aClonMat <- read.table("adeno_clonMat.txt", sep="\t", header=T, check.names=F, row.names=1)
aClonMatB <- read.table("adeno_clonMat_perBiopsy.txt", sep="\t", header=T, check.names=F, row.names=1)
cClonMat <- read.table("carcino_clonMat.txt", sep="\t", header=T, check.names=F, row.names=1)
cClonMatB <- read.table("carcino_clonMat_perBiopsy.txt", sep="\t", header=T, check.names=F, row.names=1)

library(parallel)
s="COAD"
intnf = bestCombs[s,1]
epinf = bestCombs[s,2]
maxland = 50
medFit <- unlist(lapply(bestFitL, median, na.rm=T))
names(medFit) <- names(bestFitL)
genes <- names(socRatL[[s]])[1:min(length(socRatL[[s]]), maxland)]
fMat <- read.table(paste("../tcga/epiFactMat_", s, ".txt", sep=""),
                   sep="\t", header=T, row.names=1, check.names=F)
fMat <- fMat[genes,]
sMat <- read.table(paste("../tcga/clonMat_", s, ".txt", sep=""),
                   sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
nmut <- apply(sMat, 1, sum)
selectMut <- c(selectMutL[[s]], cnSelectL[[s]][rownames(sMat)[grep("_", rownames(sMat))]])

## Adeno
agenes <- intersect(genes, rownames(aClonMat))
socRat <- (1 / socRatL[[s]][agenes]) ^ socPowL[[s]][agenes]
names(socRat) <- agenes
ssMat <- aClonMat[agenes,]
mutPerSamp <- apply(ssMat, 2, sum)
geneComb <- apply(ssMat, 2, function(x){paste(rownames(ssMat)[which(x==1)], collapse=":")})
adeFit <- unlist(mclapply(geneComb, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, sf=intnf, ef=epinf, mc.cores=8))
adeFitN <- adeFit / median(medFit[s], na.rm=T)

ssMatB <- aClonMatB[agenes,]
mutPerSampB <- apply(ssMatB, 2, sum)
geneCombB <- apply(ssMatB, 2, function(x){paste(rownames(ssMat)[which(x==1)], collapse=":")})
adeFitB <- unlist(mclapply(geneCombB, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, sf=intnf, ef=epinf, mc.cores=8))
adeFitNB <- adeFitB / median(medFit[s], na.rm=T)

## Carcino
cgenes <- intersect(genes, rownames(cClonMat))
socRat <- (1 / socRatL[[s]][cgenes]) ^ socPowL[[s]][cgenes]
names(socRat) <- cgenes
cssMat <- cClonMat[cgenes,]
cmutPerSamp <- apply(cssMat, 2, sum)
cgeneComb <- apply(cssMat, 2, function(x){paste(rownames(cssMat)[which(x==1)], collapse=":")})
carFit <- unlist(mclapply(cgeneComb, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, sf=intnf, ef=epinf, mc.cores=8))
carFitN <- carFit / median(medFit[s], na.rm=T)

cssMatB <- cClonMatB[cgenes,]
cmutPerSampB <- apply(cssMatB, 2, sum)
cgeneCombB <- apply(cssMatB, 2, function(x){paste(rownames(cssMat)[which(x==1)], collapse=":")})
carFitB <- unlist(mclapply(cgeneCombB, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, sf=intnf, ef=epinf, mc.cores=8))
carFitNB <- carFitB / median(medFit[s], na.rm=T)

rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")
pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- names(socPowL)
colSet <- rainbow10[1:length(socPowL)]
names(colSet) <- names(socPowL)
library(beeswarm)
pdf(paste("plots/adeno_fitness_",model,".pdf",sep=""))
wt <- wilcox.test(carFitN, adeFitN)
boxplot(carFitN, adeFitN, sub=sprintf("p=%.1e  (Wilcoxon test)", wt$p.value), main="Clonal mutations per adenoma",
          names=c("Carcinomas", "Adenomas"), ylab="Normalised fitness", outline=F,
          ylim=c(0, max(c(carFitN, adeFitN), na.rm=T)))
beeswarm(c(carFitN, adeFitN) ~ rep(1:2, c(length(carFitN), length(adeFitN))), add=T, cex=1.5,
         pch=pchSet[s], pwcol=rep(rainbow[c(1,7)], c(length(carFitN), length(adeFitN))), corral="wrap")

wt <- wilcox.test(carFitNB, adeFitNB)
boxplot(carFitNB, adeFitNB, sub=sprintf("p=%.1e  (Wilcoxon test)", wt$p.value), main="Clonal mutations per biopsy",
          names=c("TCGA COAD", "Adenomas"), ylab="Normalised fitness", outline=F,
          ylim=c(0, max(c(carFitNB, adeFitNB), na.rm=T)))
beeswarm(c(carFitNB, adeFitNB) ~ rep(1:2, c(length(carFitNB), length(adeFitNB))), add=T, cex=1.5,
         pch=pchSet[s], pwcol=rep(rainbow[c(1,7)], c(length(carFitNB), length(adeFitNB))), corral="wrap")
dev.off()

## Fitness in other sets
sets <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "SKCM")
s="COAD"
osets <- setdiff(sets, s)
fitMat <- matrix(nr=length(adeFitN), nc=length(osets))
fitMatB <- matrix(nr=length(adeFitNB), nc=length(osets))
rownames(fitMat) <- names(adeFitN)
rownames(fitMatB) <- names(adeFitNB)
colnames(fitMat) <- colnames(fitMatB) <- c(osets)
library(parallel)
for (s in osets) {
  print(s)
  ssf <- bestCombs[s,1]
  sef <- bestCombs[s,2]
  genes <- topL[[s]][1:min(length(topL[[s]]), maxland)]
  fMat <- read.table(paste("../tcga/epiFactMat_", s, ".txt", sep=""),
                     sep="\t", header=T, row.names=1, check.names=F)
  fMat <- fMat[genes,]
  sMat <- read.table(paste("../tcga/clonMat_", s, ".txt", sep=""),
                     sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
  nmut <- apply(sMat, 1, sum)
  selectMut <- c(selectMutL[[s]], cnSelectL[[s]][rownames(sMat)[grep("_", rownames(sMat))]])
  genes <- intersect(genes, rownames(aClonMat))
  socRat <- (1 / socRatL[[s]][genes]) ^ socPowL[[s]][genes]
  names(socRat) <- genes

  ssMat <- aClonMat[genes,]
  mutPerSamp <- apply(ssMat, 2, sum)
  geneComb <- apply(ssMat, 2, function(x){paste(rownames(ssMat)[which(x==1)], collapse=":")})
  sFit <- unlist(mclapply(geneComb, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, sf=ssf, ef=sef, mc.cores=8))
  sFitN <- sFit / median(medFit[s], na.rm=T)
  fitMat[,s] <- sFitN

  ssMatB <- aClonMatB[genes,]
  mutPerSampB <- apply(ssMatB, 2, sum)
  geneCombB <- apply(ssMatB, 2, function(x){paste(rownames(ssMat)[which(x==1)], collapse=":")})
  sFitB <- unlist(mclapply(geneCombB, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, sf=ssf, ef=sef, mc.cores=8))
  sFitNB <- sFitB / median(medFit[s], na.rm=T)
  fitMatB[,s] <- sFitNB
}

meanOther <- apply(fitMat, 1, mean)
meanOtherB <- apply(fitMatB, 1, mean)

s="COAD"
pdf(paste("plots/adeno_other_fitness_",model,".pdf",sep=""))
par(mar=c(6,4,1,1))
ptt <- t.test(adeFitN, meanOther, paired=T)
boxplot(cbind(adeFitN, meanOther), outline=F, ylim=c(0, max(c(adeFitN, meanOther))), names=c("COAD", "Other sets (mean)"),
        ylab="Normalised Fitness", sub=sprintf("p=%.1e (paired t-test)", ptt$p.value))
beeswarm(c(adeFitN, meanOther) ~ rep(1:2, rep(length(sFitN), 2)), pch=pchSet[s],
         pwcol=rep(rainbow[c(1,7)], c(length(adeFitN), length(meanOther))), add=T, cex=2)

ptt <- t.test(adeFitNB, meanOtherB, paired=T)
boxplot(cbind(adeFitNB, meanOtherB), outline=F, ylim=c(0, max(c(adeFitNB, meanOtherB))), names=c("COAD", "Other sets (mean)"),
        ylab="Normalised Fitness", sub=sprintf("p=%.1e (paired t-test)", ptt$p.value))
beeswarm(c(adeFitNB, meanOtherB) ~ rep(1:2, rep(length(sFitNB), 2)), pch=pchSet[s],
         pwcol=rep(rainbow[c(1,7)], c(length(adeFitNB), length(meanOtherB))), add=T, cex=2)
dev.off()
