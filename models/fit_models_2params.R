## Script 2
## Try all combinations of factors (listed below) and determine the best combination for each tumour type.
## Run twice, once with the combined_prod model (default) and once with the combined_mean model.
## Careful to mclapply calls, change mc.cores parameter if not enough processors are available.
##
## Notes:
## socnf= Sociability numeric factor: weight for self-sufficiency * self
## epinf= Epistatic numeric factor: weight for interactions

rm(list=ls())
gc()
model="combined_prod"
##model="combined_mean"

## multiplicative factors assayed to weight components
multf <- c(1/c(seq(100, 10, -10),9:1), 2:9, seq(10, 100, 10))
## Colours
rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")

## Data on select CN events
load("../misc/cnSelectL.RData")
## list of selection values 
load("../tcga/selectMutL.RData")

#### Objective function ####
## Normalise all scores by the median, calculate MAD (weighted average over all sets)

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

getFitStr <- function(panelStr, sel, socRat, ints, sf, ef) {
  if(panelStr=="")
    return(0)
  subpanel=unlist(strsplit(panelStr, ":"))
  return(sum(unlist(lapply(subpanel, function(g){(sel[g] * socRat[g]) * sf + getEpi(ints, g, subpanel) * ef}))))
}

## Most clonally mutated genes
load("../tcga/topL.RData")
## list of ratio of mean observed n partners compared to mean expected n partners
load("../tcga/socRatL.RData")
## list of multer values for the observed ratios
load("../tcga/socPowL.RData")

## Loop all sets
library(parallel)
## Maximum number of drivers in landscape
maxland = 50
allFitL <- list()
facTab <- matrix(NA, nc=2, nr=length(multf)*length(multf))
madTab <- matrix(NA, nc=length(socPowL), nr=length(multf)*length(multf))
sdTab <- matrix(NA, nc=length(socPowL), nr=length(multf)*length(multf))
colnames(facTab) <- c("IntFactor", "EpiFactor")
colnames(madTab) <- colnames(sdTab) <- c(names(socPowL))
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
  fMat <- fMat[genes,genes]
  geneComb <- apply(clonMat, 2, function(x){paste(rownames(clonMat)[which(x==1)], collapse=":")})

  allFitL[[s]] <- matrix(nc=ncol(clonMat), nr=nrow(madTab))
  colnames(allFitL[[s]]) <- colnames(clonMat)
  i=0
  for (socnf in (multf)) {
    for (epinf in (multf)) {
      i=i+1
      if (s == names(socPowL)[1])
        facTab[i,] <- c(socnf, epinf)
      allFit <- unlist(mclapply(geneComb, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, sf=socnf, ef=epinf, mc.cores=8))
      m <- median(allFit, na.rm=T)
      if (m == 0)
        m <- min(allFitN[which(!is.na(allFitN & allFitN > 0))], na.rm=T)
      allFitN <- allFit / m
      madTab[i,s] <- mean(abs(allFitN - 1), na.rm=T)
      sdTab[i,s] <- sd(allFitN, na.rm=T)
      allFitL[[s]][i,] <- allFit
    }
  }
}

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

bestCombs <- facTab[bestIdx2,,drop=F]
rownames(bestCombs) <- names(allFitL)
bestCombs
save(bestCombs, file=paste("bestCombs_",model,".RData",sep=""))

pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- names(socPowL)
colSet <- rainbow10[1:length(socPowL)]
names(colSet) <- names(socPowL)

library(gplots)
library(beeswarm)
pdf(paste("plots/factor_heatmaps_",model,".pdf",sep=""))
for (s in names(allFitL)) {
  ## madMat <- matrix(madTabL[[s]], nc=length(multf), nr=length(multf), byrow=F)
  ## colnames(madMat) <- sprintf("%.2f", multf)
  ## rownames(madMat) <- sprintf("%.2f", multf)
  ## heatmap.2(madMat, trace="none", dendrogram="none", Rowv=NA, Colv=NA, main=paste(s, "MAD"), symkey=F, col=bluered(200),
  ##           ylab="Epistatic Factor", xlab="Intrinsic Factor")
  sdMat <- matrix(sdTabL[[s]], nc=length(multf), nr=length(multf), byrow=F)
  colnames(sdMat) <- sprintf("%.2f", multf)
  rownames(sdMat) <- sprintf("%.2f", multf)
  heatmap.2(sdMat, trace="none", dendrogram="none", Rowv=NA, Colv=NA, main=paste(s, "Std dev"), symkey=F, col=bluered(200),
            ylab="Epistatic Factor", xlab="Intrinsic Factor")  
}
dev.off()

load(paste("bestCombs_",model,".RData", sep=""))
pdf(paste("plots/factor_inference_",model,".pdf",sep=""))
## plot(log10(bestCombs), pch=pchSet[names(allFitL)], col=colSet[names(allFitL)], xlab="Optimal Epistatic Factor (log2)",
##     ylab="Optimal Intrinsic Factor (log10)", cex=2, main="MAD", xlim=log10(c(min(multf), max(multf))),
##     ylim=log10(c(min(multf), max(multf))))
## legend("topright", legend=names(allFitL), col=colSet[names(allFitL)], pch=pchSet[names(allFitL)])
plot(log10(bestCombs), pch=pchSet[names(allFitL)], col=colSet[names(allFitL)], xlab="Optimal Epistatic Factor (log2)",
     ylab="Optimal Intrinsic Factor (log10)", cex=2, main="Standard deviation", xlim=log10(c(min(multf), max(multf))),
     ylim=log10(c(min(multf), max(multf))))
legend("topright", legend=names(allFitL), col=colSet[names(allFitL)], pch=pchSet[names(allFitL)])
dev.off()

## Calculate share of both indices
getFitComponents <- function(panelStr, sel, socRat, ints, sf, ef) {
  subpanel=unlist(strsplit(panelStr, ":"))
  ##print(subpanel)
  if (length(subpanel) == 0)
    return(NA)
  ##return(sum(unlist(lapply(subpanel, function(g){occur[g] * (socRat[g] * sf + getEpi(ints, g, subpanel) * ef)}))))
  return(unlist(lapply(subpanel, function(g){c((sel[g] * socRat[g]) * sf, getEpi(ints, g, subpanel) * ef)})))
}

library(parallel)
library(gplots)
library(beeswarm)
library(coin)
nDrivL <- bestFitL <- bestFitNL <- bestIntL <- bestEpiL <- list()
pchmut <- c(1,0,2,5,16,15,17,18)
if (length(pchmut) < maxland) {
  pchmut <- c(pchmut, rep(18, maxland - length(pchmut)))
}
names(pchmut) <- as.character(1:length(pchmut))
pchmut["0"] <- 3
allMat <- matrix(nr=0, ncol=7)
colnames(allMat) <- c("Set", "IntScore", "EpiScore", "AllFit", "AllFitN", "IntPerc", "MutPerSamp")
spercL <- list()

## Calculate fitness using best weigths and save it.
## Calculate correlations between components, number of drivers and the
## share of each component in the fitness score of each sample.
pdf(paste("plots/component_prevalence_", model, ".pdf", sep=""))
for (s in names(socPowL)) {
  socnf = bestCombs[s,"IntFactor"]
  epinf = bestCombs[s,"EpiFactor"]
  print(s)
  genes <- topL[[s]][1:min(length(topL[[s]]), maxland)]
  clonMat <- read.table(paste("../tcga/clonMat_", s, ".txt", sep=""), sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
  nmut <- apply(clonMat, 1, sum)
  selectMut <- c(selectMutL[[s]], cnSelectL[[s]][rownames(clonMat)[grep("_", rownames(clonMat))]])
  clonMat <- clonMat[genes,]
  mutPerSamp <- apply(clonMat, 2, function(x){length(which(x>0))})
  socRat <- (1 / socRatL[[s]][genes]) ^ socPowL[[s]][genes]
  names(socRat) <- genes
  fMat <- read.table(paste("../tcga/epiFactMat_", s, ".txt", sep=""), sep="\t", header=T, row.names=1, check.names=F)
  fMat <- fMat[genes,]
  geneComb <- apply(clonMat, 2, function(x){paste(rownames(clonMat)[which(x==1)], collapse=":")})
  pchmut[as.character(unique(mutPerSamp[which(mutPerSamp > maxland)]))] <- pchmut[as.character(maxland)]
  nDrivL[[s]] <- mutPerSamp
  bestFitL[[s]] <- allFit <- unlist(mclapply(geneComb, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, sf=socnf, ef=epinf, mc.cores=8))
  bestFitNL[[s]] <- allFitN <- allFit / median(allFit, na.rm=T)
  geneL <- lapply(geneComb, function(x){unlist(strsplit(x, ":"))})
  fitcomp <- lapply(geneComb, getFitComponents, sel=selectMut, socRat=socRat, ints=fMat, sf=socnf, ef=epinf)
  spercL[[s]] <- unlist(lapply(fitcomp, function(x){
    if (length(x) < 2)
      return(NA)
    s=x[seq(1, (length(x)-1), 2)]
    sum(s)/sum(x)
  }))

  both <- lapply(fitcomp, function(x){
    if (length(x) < 2)
      return(c(NA,NA))
    s=x[seq(1, length(x)-1, 2)]
    e=x[seq(2, length(x), 2)]
    return(c(sum(s),sum(e)))
  })
  bothMat <- matrix(unlist(both), nc=2, byrow = T)
  bestIntL[[s]] <- bothMat[,1]
  bestEpiL[[s]] <- bothMat[,2]
  setMat <- cbind(bothMat, allFit, allFitN, spercL[[s]], mutPerSamp)
  setMat <- setMat[which(!is.na(allFit)),]
  allMat <- rbind(allMat, cbind(rep(s, nrow(setMat)), setMat))

  it <- independence_test(allFitN ~ mutPerSamp)
  boxplot(allFitN ~ mutPerSamp, ylab="Normalised Fitness", xlab="Number of clonal drivers", outline=F, ylim=c(0, max(allFitN, na.rm=T)),
          sub=sprintf("p=%.1e (independence test)", pvalue(it)))
  beeswarm(allFitN ~ mutPerSamp, pch=pchSet[s], col=colSet[s], add=T, corral="wrap")
}

boxplot(spercL, ylab="Percentage of Intrinsic component in Fitness score", outline=F,
        ylim=c(min(unlist(spercL), na.rm=T), max(unlist(spercL), na.rm=T)), las=2)
beeswarm(unlist(spercL) ~ rep(1:length(spercL), unlist(lapply(spercL, length))),
         add=T, corral="wrap", cex=0.5,
         pwpch=pchSet[rep(names(spercL), unlist(lapply(spercL, length)))],
         pwcol=colSet[rep(names(spercL), unlist(lapply(spercL, length)))])

allFit <- unlist(bestFitNL)
allNDriv <- unlist(nDrivL)
allPch <- pchSet[rep(names(socPowL), unlist(lapply(bestFitNL, length)))]
allCol <- colSet[rep(names(socPowL), unlist(lapply(bestFitNL, length)))]
ct <- cor.test(allInt, allFit)
it <- independence_test(allFit ~ allNDriv)
boxplot(allFit ~ allNDriv, ylab="Normalised Fitness", xlab="Number of clonal drivers", outline=F, ylim=c(0, max(allFit, na.rm=T)),
        sub=sprintf("p=%.1e (independence test)", pvalue(it)), las=2)
beeswarm(allFit ~ allNDriv, pwpch=allPch, pwcol=allCol, add=T, corral="wrap", cex=0.5)
legend("topright", legend=names(bestFitNL), col=colSet[names(bestFitNL)], pch=pchSet[names(bestFitNL)])
dev.off()

save(bestFitNL, file=paste("bestFitNL_",model,".RData",sep=""))
save(bestFitL, file=paste("bestFitL_",model,".RData",sep=""))
