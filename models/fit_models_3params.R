## Script 3
## Try all combinations of factors (listed below) and determine the best combination for each tumour type.
## Run twice, once with the separate_prod model (default) and once with the separate_mean model.
## Careful to mclapply calls, change mc.cores parameter if not enough processors are available.
##
## Notes:
## Occnf= Occurrence numeric factor: weight for selection strength
## Socnf= Sociability numeric factor: weight for self-sufficiency
## Epinf= Epistatic numeric factor: weight for interactions

rm(list=ls())
gc()
model="separate_prod"
##model="separate_mean"

#### Parameters
multf <- c(1/c(seq(100, 10, -10),9:1), 2:9, seq(10, 100, 10))

## Data on select CN events
load("../tcga/cnSelectL.RData")
## list of selection values 
load("../tcga/selectMutL.RData")

#### Objective function ####
## Normalise all scores by the median, calculate MAD.
## 1) average over all sets. 2) average of averages per set.

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

getFitStr <- function(panelStr, sel, socRat, ints, of, sf, ef) {
  subpanel=unlist(strsplit(panelStr, ":"))
  ##print(subpanel)
  return(sum(unlist(lapply(subpanel, function(g){(sel[g] * of) + (socRat[g] * sf) + (getEpi(ints, g, subpanel) * ef)}))))
}

## Most clonally mutated genes
load("../tcga/topL.RData")
## list of ratio of mean observed n partners compared to mean expected n partners
load("../tcga/socRatL.RData")
## list of power values for the observed ratios
load("../tcga/socPowL.RData")

library(parallel)
## Loop all sets
maxland = 50
allFitL <- list()
facTab <- matrix(NA, nc=3, nr=length(multf)*length(multf)*length(multf))
madTab <- matrix(NA, nc=length(socPowL), nr=length(multf)*length(multf)*length(multf))
colnames(facTab) <- c("OccFactor", "SocFactor", "EpiFactor")
colnames(madTab) <- colnames(madTab) <- c(names(socPowL))
for (s in names(socPowL)) {
  print(s)
  genes <- topL[[s]][1:min(length(topL[[s]]), maxland)]
  clonMat <- read.table(paste("../tcga/clonMat_", s, ".txt", sep=""), sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
  nmut <- apply(clonMat, 1, sum)
  selectMut <- c(selectMutL[[s]], cnSelectL[[s]][rownames(clonMat)[grep("_", rownames(clonMat))]])
  clonMat <- clonMat[genes,]
  mutPerSamp <- apply(clonMat, 2, sum)
  socRat <- (1 / socRatL[[s]][genes]) ^ socPowL[[s]][genes]
  names(socRat) <- genes
  fMat <- read.table(paste("../tcga/epiFactMat_", s, ".txt", sep=""), sep="\t", header=T, row.names=1, check.names=F)
  fMat <- fMat[genes,]
  geneComb <- apply(clonMat, 2, function(x){paste(rownames(clonMat)[which(x==1)], collapse=":")})

  allFitL[[s]] <- matrix(nc=ncol(clonMat), nr=nrow(madTab))
  colnames(allFitL[[s]]) <- colnames(clonMat)
  i=0
  for (occnf in multf) {
    print(occnf)
    for (socnf in multf) {
      print(paste("", socnf))
      for (epinf in multf) {
        ##print(paste(" ", epinf))
        i=i+1
        if (s == names(socPowL)[1])
          facTab[i,] <- c(occnf, socnf, epinf)
        allFit <- unlist(mclapply(geneComb, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, of=occnf, sf=socnf, ef=epinf, mc.cores=8))
        allFitN <- allFit / median(allFit, na.rm=T)
        madTab[i,s] <- mean(abs(allFitN - 1), na.rm=T)
        allFitL[[s]][i,] <- allFit
      }
    }
  }
}
save(allFitL, file=paste("allFitL_", model, ".RData", sep=""))

bestIdx <- bestIdx2 <- c()
madTabL <- sdTabL <- list()
for (s in names(allFitL)) {
  madTab <- sdTab <- c()
  for (i in 1:nrow(facTab)) {
    m <- median(allFitL[[s]][i,], na.rm=T)
    if (m == 0)
      m <- mean(allFitL[[s]][i,])      
    tmp <- allFitL[[s]][i,] / m
    madTab[i] <- mean(abs(tmp-1), na.rm=T)
    sdTab[i] <- sd(tmp-1, na.rm=T)
  }
  madTabL[[s]] <- madTab
  sdTabL[[s]] <- sdTab
  bestIdx[s] <- which(madTab == min(madTab))[1]
  bestIdx2[s] <- which(sdTab == min(sdTab))[1]
}

## Best Occur, Soc, Epi factors
bestCombs <- facTab[bestIdx2,]
rownames(bestCombs) <- names(allFitL)
save(bestCombs, file=paste("bestCombs_",model,".RData",sep=""))

### Plot for best values each time
pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- names(socPowL)
colSet <- rainbow10[1:length(socPowL)]
names(colSet) <- names(socPowL)

library(gplots)
library(beeswarm)
pdf(paste("plots/factor_heatmaps_",model,".pdf",sep=""))
for (s in rownames(bestCombs)) {
  besto <- bestCombs[s,"OccFactor"]
  bests <- bestCombs[s,"SocFactor"]
  beste <- bestCombs[s,"EpiFactor"]

  ## Best occurrence factor
  sdMat <- matrix(sdTabL[[s]][which(facTab[,"OccFactor"] == besto)], nc=length(multf), nr=length(multf), byrow=F)
  colnames(sdMat) <- rownames(sdMat) <- sprintf("%.2f", multf)
  heatmap.2(sdMat, trace="none", dendrogram="none", Rowv=NA, Colv=NA, main=sprintf("%s SelFactor=%.2f", s, besto),
            symkey=F, col=bluered(200), ylab="Epistatic Factor", xlab="Sociable Factor")

  sdMat <- matrix(sdTabL[[s]][which(facTab[,"SocFactor"] == bests)], nc=length(multf), nr=length(multf), byrow=F)
  colnames(sdMat) <- rownames(sdMat) <- sprintf("%.2f", multf)
  heatmap.2(sdMat, trace="none", dendrogram="none", Rowv=NA, Colv=NA, main=sprintf("%s SocFactor=%.2f", s, bests),
            symkey=F, col=bluered(200), ylab="Epistatic Factor", xlab="Selection Factor")

  sdMat <- matrix(sdTabL[[s]][which(facTab[,"EpiFactor"] == beste)], nc=length(multf), nr=length(multf), byrow=F)
  colnames(sdMat) <- rownames(sdMat) <- sprintf("%.2f", multf)
  heatmap.2(sdMat, trace="none", dendrogram="none", Rowv=NA, Colv=NA, main=sprintf("%s EpiFactor=%.2f", s, beste),
            symkey=F, col=bluered(200), ylab="Sociable Factor", xlab="Selection Factor")
}
dev.off()

library(plot3D)
load(paste("bestCombs_",model,".RData", sep=""))
pdf(paste("plots/factor_inference_",model,".pdf",sep=""))
scatter3D(x=bestCombs[,1], y=bestCombs[,2], z=bestCombs[,3], pch=pchSet[rownames(bestCombs)], colvar=1:9, col=colSet, cex=3,
          xlab="Selection", ylab="Self-sufficiency", zlab="Epistasis", colkey=F, theta=60, phi=30, xlim=c(0, max(multf)),
           ticktype= "detailed", bty = "g", alpha = 0.5, ylim=c(0, max(multf)), zlim=c(0, max(multf)))
legend("bottomright", legend=rownames(bestCombs), col=colSet[rownames(bestCombs)], pch=pchSet[rownames(bestCombs)], bty="n")
scatter3D(x=log10(bestCombs[,1]), y=log10(bestCombs[,2]), z=log10(bestCombs[,3]), ticktype= "detailed",
          bty = "g", alpha = 0.5, pch=pchSet[rownames(bestCombs)], colvar=1:9, col=colSet, cex=3,
          xlab="Selection", ylab="Self-sufficiency", zlab="Epistasis", colkey=F, theta=60, phi=30,
          xlim=log10(c(min(multf), max(multf))), ylim=log10(c(min(multf), max(multf))), zlim=log10(c(min(multf), max(multf))))
legend("bottomright", legend=rownames(bestCombs), col=colSet[rownames(bestCombs)], pch=pchSet[rownames(bestCombs)], bty="n")
scatter3D(x=bestCombs[,1], y=bestCombs[,2], z=bestCombs[,3], pch=pchSet[rownames(bestCombs)], colvar=1:9, col=colSet, ticktype= "detailed",
          bty = "g", alpha = 0.5, cex=3, xlab="Selection", ylab="Self-sufficiency", zlab="Epistasis", colkey=F,
          theta=15, phi=0, xlim=c(0, max(multf)), ylim=c(0, max(multf)), zlim=c(0, max(multf)))
legend("topright", legend=rownames(bestCombs), col=colSet[rownames(bestCombs)], pch=pchSet[rownames(bestCombs)], bty="n")
scatter3D(x=log10(bestCombs[,1]), y=log10(bestCombs[,2]), z=log10(bestCombs[,3]), ticktype= "detailed",
          bty = "g", alpha = 0.5, pch=pchSet[rownames(bestCombs)], colvar=1:9, col=colSet, cex=3,
          xlab="Selection", ylab="Self-sufficiency", zlab="Epistasis", colkey=F, theta=15, phi=0,
          xlim=log10(c(min(multf), max(multf))), ylim=log10(c(min(multf), max(multf))), zlim=log10(c(min(multf), max(multf))))
legend("topright", legend=rownames(bestCombs), col=colSet[rownames(bestCombs)], pch=pchSet[rownames(bestCombs)], bty="n")
dev.off()

## Calculate share of both indices

getFitComponents <- function(panelStr, ints, sel, socRat, of, sf, ef) {
  subpanel=unlist(strsplit(panelStr, ":"))
  ##print(subpanel)
  if (length(subpanel) == 0)
    return(NA)
  return(unlist(lapply(subpanel, function(g){c(sel[g] * of, socRat[g] * sf, getEpi(ints, g, subpanel) * ef)})))
}

library(coin)
rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")
pchmut <- c(1,0,2,5,16,15,17,18)
if (length(pchmut) < maxland) {
  pchmut <- c(pchmut, rep(18, maxland - length(pchmut)))
}
pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- names(socPowL)
colSet <- rainbow10[1:length(socPowL)]
names(colSet) <- names(socPowL)
names(pchmut) <- as.character(1:length(pchmut))
pchmut["0"] <- 3

allMat <- matrix(nr=0, ncol=11)
colnames(allMat) <- c("Set", "OccScore", "SocScore", "EpiScore", "AllFit", "AllFitN", "OccPerc", "SocPerc", "EpiPerc", "MutPerSamp", "FitCol")
bestFitL <- bestFitNL <- nDrivL <- list()
for (s in names(socPowL)) {
  print(s)
  occnf <- bestCombs[s,"OccFactor"]
  socnf <- bestCombs[s,"SocFactor"]
  epinf <- bestCombs[s,"EpiFactor"]
  genes <- topL[[s]][1:min(length(topL[[s]]), maxland)]
  clonMat <- read.table(paste("../clonMat_", s, ".txt", sep=""), sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
  clonMat <- clonMat[genes,]
  nmut <- apply(clonMat, 1, sum)
  mutPerSamp <- apply(clonMat, 2, sum)
  nDrivL[[s]] <- mutPerSamp
  selectMut <- c(selectMutL[[s]], cnSelectL[[s]][rownames(clonMat)[grep("_", rownames(clonMat))]])
  socRat <- (1 / socRatL[[s]][genes]) ^ socPowL[[s]][genes]
  names(socRat) <- genes
  fMat <- read.table(paste("../epiFactMat_", s, ".txt", sep=""), sep="\t", header=T, row.names=1, check.names=F)
  fMat <- fMat[genes,]

  geneComb <- apply(clonMat, 2, function(x){paste(rownames(clonMat)[which(x==1)], collapse=":")})
  bestFitL[[s]] <- allFit <- allFitL[[s]][bestIdx2[s],]
  bestFitNL[[s]] <- allFitN <- allFit / median(allFitL[[s]][bestIdx2[s],])
}
save(bestFitNL, file=paste("bestFitNL_",model,".RData",sep=""))
save(bestFitL, file=paste("bestFitL_",model,".RData",sep=""))

pdf(paste("plots/component_prevalence_",model,".pdf", sep=""))
for (s in names(socPowL)) {
  print(s)
  occnf <- bestCombs[s,"OccFactor"]
  socnf <- bestCombs[s,"SocFactor"]
  epinf <- bestCombs[s,"EpiFactor"]
  genes <- topL[[s]][1:min(length(topL[[s]]), maxland)]
  clonMat <- read.table(paste("../clonMat_", s, ".txt", sep=""), sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
  clonMat <- clonMat[genes,]
  geneComb <- apply(clonMat, 2, function(x){paste(rownames(clonMat)[which(x==1)], collapse=":")})
  fitcomp <- lapply(geneComb, getFitComponents, ints=fMat, sel=selectMut, socRat=socRat, of=occnf, sf=socnf, ef=epinf)
  fitperc <- matrix(unlist(lapply(fitcomp, function(x){
    if (length(x) < 3)
      return(c(NA,NA,NA))
    c(sum(x[seq(1, length(x)-2, 3)]) / sum(x),
      sum(x[seq(2, length(x)-1, 3)]) / sum(x),
      sum(x[seq(3, length(x), 3)]) / sum(x))})), nc=3, byrow=T)
  
  boxplot(fitperc, ylab="Percentage of component in Fitness score", outline=F,
          ylim=c(min(fitperc, na.rm=T), max(fitperc, na.rm=T)), main=s, names=c("Selection", "Social", "Epistatic"))
  beeswarm(as.vector(fitperc) ~ rep(1:ncol(fitperc), rep(nrow(fitperc), ncol(fitperc))), add=T, pch=pchSet[s], col=colSet[s],
           corral="wrap", cex=0.5)

  it <- independence_test(allFitN ~ mutPerSamp)
  boxplot(allFitN ~ mutPerSamp, ylab="Normalised Fitness", xlab="Number of clonal drivers", outline=F, ylim=c(0, max(allFitN, na.rm=T)),
          sub=sprintf("p=%.1e (independence test)", pvalue(it)))
  beeswarm(allFitN ~ mutPerSamp, pch=pchSet[s], col=colSet[s], add=T, corral="wrap")
}
allFit <- unlist(bestFitNL)
allNDriv <- unlist(nDrivL)
allPch <- pchSet[rep(names(socPowL), unlist(lapply(bestFitNL, length)))]
allCol <- colSet[rep(names(socPowL), unlist(lapply(bestFitNL, length)))]
allCex <- unlist(nDrivL)
allCex[allCex > 2.5] <- 2.5
it <- independence_test(allFit ~ allNDriv)
boxplot(allFit ~ allNDriv, ylab="Normalised Fitness", xlab="Number of clonal drivers", outline=F, ylim=c(0, max(allFit, na.rm=T)),
        sub=sprintf("p=%.1e (independence test)", pvalue(it)), las=2)
beeswarm(allFit ~ allNDriv, pwpch=allPch, pwcol=allCol, add=T, corral="wrap", cex=0.5)
legend("topright", legend=names(bestFitNL), col=colSet[names(bestFitNL)], pch=pchSet[names(bestFitNL)])
dev.off()
