rm(list=ls())
gc()
model="combined_prod"

#### Parameters
##intf <- c(1/(10:1), 2:10)
powf <- c(1/c(5:3, seq(2, 1, -0.2)), seq(1.2, 2, 0.2), 3:5)
##powf <- seq(.1, 2, 0.1)
multf <- c(1/c(seq(100, 10, -10),9:1), 2:9, seq(10, 100, 10))
##modif <- c("self", "log", "square")
op <- c("+", "x")
rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")

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
load("../tcga/cnSelectL.RData")

## util functions
getEpi <- function(ints, g, subpanel) {
  partners <- setdiff(subpanel,g)
  if (length(partners) > 0) {
    if (length(grep("_mean_perSet$", model)) > 0) {
      return(mean(unlist(ints[g,setdiff(subpanel,g)])))
    } else {
      return(prod(ints[g,setdiff(subpanel,g)]))
    }
  } else {
    return(1)
  }
}

getFit <- function(subpanel, sel, socRat, ints, sf, ef) {
  if(length(subpanel) == 0)
    return(0)
  return(sum(unlist(lapply(subpanel, function(g){(sel[g] * socRat[g]) * sf + getEpi(ints, g, subpanel) * ef}))))
}

getFitStr <- function(panelStr, sel, socRat, ints, sf, ef) {
  if(panelStr=="")
    return(0)
  subpanel=unlist(strsplit(panelStr, ":"))
  return(sum(unlist(lapply(subpanel, function(g){(sel[g] * socRat[g]) * sf + getEpi(ints, g, subpanel) * ef}))))
}

getExpectMut <- function(g, totMut, len) {
  (len[g] / sum(len, na.rm=T)) * totMut
}

getCombFit <- function(allComb, genes, sel, socRat, ints, sf, ef, mf) {
  combFit <- unlist(apply(allComb, 2, function(gidx){getFit(genes[gidx], sel, socRat, ints, sf, ef) / mf}))
  names(combFit) <- apply(allComb, 2, function(x){paste(genes[x],collapse=":")})
  return(combFit)
}

getParents <- function(x) {
  splt <- unlist(strsplit(x, split=":"))
  combs <- combn(length(splt),length(splt)-1)
  return(unlist(apply(combs,2, function(idx,l){paste(l[idx], collapse=":")},l=splt)))
}

getEdges <- function(n, combFit) {
  tab <- c()
  print(n)
  geneCombs <- names(combFit[[n]])
  parents <- mclapply(geneCombs, getParents, mc.cores=10)
  return(cbind(unlist(parents), rep(geneCombs, rep(length(parents[[1]]), length(geneCombs)))))
}

#################################
####     Create networks     ####
#################################

library(parallel)
maxland = 15
load(paste("../models/bestCombs_",model,".RData", sep=""))
load(paste("../models/bestFitNL_",model,".RData", sep=""))
load(paste("../models/bestFitL_",model,".RData", sep=""))
medFit <- unlist(lapply(bestFitL, median, na.rm=T))
for (s in names(bestFitL)) {
  print(s)
  intnf = bestCombs[s,1]
  epinf = bestCombs[s,2]
  genes <- topL[[s]][1:min(length(topL[[s]]), maxland)]
  socRat <- (1 / socRatL[[s]][genes]) ^ socPowL[[s]][genes]
  names(socRat) <- genes
  fMat <- read.table(paste("../tcga/epiFactMat_", s, ".txt", sep=""), sep="\t", header=T, row.names=1, check.names=F)
  fMat <- fMat[genes,]
  sMat <- read.table(paste("../tcga/clonMat_", s, ".txt", sep=""),
                     sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
  sMat <- sMat[genes,]
  nmut <- apply(sMat, 1, sum)
  selectMut <- c(selectMutL[[s]], cnSelectL[[s]][rownames(sMat)[grep("_", rownames(sMat))]])
  
  geneComb <- mclapply(1:length(genes), combn, x=length(genes), mc.cores=8)

  ## Get nodes
  combFit <- mclapply(geneComb, getCombFit, genes=genes, sel=selectMut, socRat=socRat, ints=fMat,
                      sf=intnf, ef=epinf, mf=medFit[s], mc.cores=8)
  ## combNsamp <- mclapply(geneComb, combNsample, genes=genes, presMat=presMat, gprob=sprob, mc.cores=12)
  ## save(combNsamp, file=sprintf("%s_%d_combNsamp5.RData", s, maxland))
  nodeTab <- cbind(c("Node0",names(unlist(combFit))), c(0,unlist(combFit)))
  colnames(nodeTab) <- c("Genes", "Fitness")
  nodeRank <- unlist(lapply(nodeTab[,"Genes"], function(x){length(unlist(strsplit(x,":")))}))
  nodeRank[1] <- 0
  nodeOrder <- unlist(lapply(0:max(nodeRank), function(x, tab, nrank){order(tab[which(nrank == x),"Genes"])}, nodeTab, nodeRank))
  nodeTab <- cbind(nodeTab, nodeRank, nodeOrder)
  colnames(nodeTab) <- c("Genes", "Fitness", "Rank", "Order")
  write.table(nodeTab, file=sprintf("%s_%d_combFit_nodes.txt", s, maxland), sep="\t", row.names=F, quote=F)

  ## Get edges
  edgeTab <- matrix(cbind(rep("Node0", length(genes)), genes, nmut/ncol(sMat)), nc=3)
  edgeL <- lapply(2:length(combFit), getEdges, combFit=combFit)
  for (i in 1:length(edgeL)) {
    tmp <- matrix(as.vector(edgeL[[i]]), nc=2)
    edgeTab <- rbind(edgeTab, cbind(tmp, apply(tmp, 1, function(x){
      l1 <- unlist(strsplit(x[1], split=":"))
      l2 <- unlist(strsplit(x[2], split=":"))
      nmut[setdiff(l2,l1)] / ncol(sMat)
    })))
  }
  colnames(edgeTab) <- c("Parent", "Offspring", "Weight")
  write.table(edgeTab, file=sprintf("%s_%d_combFit_edges.txt", s, maxland), sep="\t", row.names=F, quote=F)
}

#################################
####    Analyse networks    #####
#################################

getSubFitness <- function(r, nodeTab, edgeTab, parents=NULL, maxr=NA) {
  if (is.na(maxr))
    maxr = max(nodeTab$Rank)
  if (is.null(parents)) {
    rnodes <- nodeTab$Genes[which(nodeTab$Rank==r)]
  } else {
    rnodes <- edgeTab$Offspring[which(edgeTab$Parent %in% parents)]
  }
  if (r == maxr) {
    return(nodeTab$Fitness[which(nodeTab$Genes %in% rnodes)])
  } else {
    return(c(nodeTab$Fitness[which(nodeTab$Genes %in% rnodes)], getSubFitness(r+1, nodeTab, edgeTab, rnodes)))
  }
}

library(parallel)
panelSize <- maxland
sets <- names(bestFitL)
maxFitL <- rankFitL <- maxRankFitL <- rankImproveL <- rankImprovePL <- list()
print(panelSize)
ps=as.character(panelSize)
sets <- gsub("^([A-Z]+)_.+$", "\\1", dir(pattern=sprintf("_%d_combFit_nodes.txt",panelSize)))
maxFitL <- rankFitL <- maxRankFitL <- rankImproveL <- rankImprovePL <- improveEffL <- list()  

for (s in sets) {
  print(s)
  ## Load net
  nodeTab <- read.table(sprintf("%s_%d_combFit_nodes.txt",s,panelSize), sep="\t", header=T, stringsAsFactors=F)
  edgeTab <- read.table(sprintf("%s_%d_combFit_edges.txt",s,panelSize), sep="\t", header=T, stringsAsFactors=F)
  maxFitL[[s]] <- max(nodeTab$Fitness)
  rankFitL[[s]] <- unlist(lapply(1:max(nodeTab$Rank), function(r){mean(nodeTab$Fitness[which(nodeTab$Rank==r)])}))
  maxRankFitL[[s]] <- unlist(lapply(1:max(nodeTab$Rank), function(r){max(nodeTab$Fitness[which(nodeTab$Rank==r)])}))
}
save(rankFitL, file=sprintf("rankFitL_%dgenes.RData", panelSize))
save(maxRankFitL, file=sprintf("maxRankFitL_%dgenes.RData", panelSize))
save(maxFitL, file=sprintf("maxFitL_%dgenes.RData", panelSize))

rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")
pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- names(socPowL)
colSet <- rainbow10[1:length(socPowL)]
names(colSet) <- names(socPowL)

library(beeswarm)
pdf(sprintf("plots/network_properties_%s_%dgenes.pdf", model, panelSize))
## Mean fitness per N drivers
plot(rankFitL[[1]], xlim=c(1,panelSize), ylim=c(min(unlist(rankFitL)),max(unlist(rankFitL))), col=colSet[1], type="l",
     xlab="Rank", ylab="Average fitness per number of drivers")
points(rankFitL[[1]], pch=pchSet[1], col=colSet[1])
abline(h=0, lty=2, col="darkgray", lwd=3)
lapply(2:length(sets),
       function(x){points(unlist(rankFitL[[x]]), pch=pchSet[x], col=colSet[x]);
         points(unlist(rankFitL[[x]]), col=colSet[x], type="l");})
legend("topleft", pch=pchSet, col=colSet, legend=sets, bty="n", cex=0.75)
plot(rankFitL[[1]] / unlist(maxFitL)[1], xlim=c(1,panelSize), ylim=c(0,1), col=colSet[1], type="l",
     xlab="Rank", ylab="Average set-normalised fitness per number of drivers")
points(rankFitL[[1]]/unlist(maxFitL)[1], pch=pchSet[1], col=colSet[1])
abline(h=0, lty=2, col="darkgray", lwd=3)
lapply(2:length(sets),
       function(x){points(unlist(rankFitL[[x]])/unlist(maxFitL)[x], pch=pchSet[x], col=colSet[x]);
         points(unlist(rankFitL[[x]])/unlist(maxFitL)[x], col=colSet[x], type="l");})
legend("topleft", pch=pchSet, col=colSet, legend=sets, bty="n", cex=0.75)
## Maximum fitness per N drivers
plot(maxRankFitL[[1]], xlim=c(1,panelSize), ylim=c(min(unlist(rankFitL)),max(unlist(maxFitL))), col=colSet[1],
     type="l", xlab="Rank", ylab="Maximum fitness per number of drivers")
points(maxRankFitL[[1]], pch=pchSet[1], col=colSet[1])
abline(h=0, lty=2, col="darkgray", lwd=3)
lapply(2:length(sets),
       function(x){points(unlist(maxRankFitL[[x]]), pch=pchSet[x], col=colSet[x]);
         points(unlist(maxRankFitL[[x]]), col=colSet[x], type="l");})
plot(maxRankFitL[[1]] / unlist(maxFitL)[1], xlim=c(1,panelSize), ylim=c(0,1), col=colSet[1], type="l",
     xlab="Rank", ylab="Maximum set-normalised fitness per number of drivers")
points(maxRankFitL[[1]]/unlist(maxFitL)[1], pch=pchSet[1], col=colSet[1])
abline(h=0, lty=2, col="darkgray", lwd=3)
lapply(2:length(sets),
       function(x){points(unlist(maxRankFitL[[x]])/unlist(maxFitL)[x], pch=pchSet[x], col=colSet[x]);
         points(unlist(maxRankFitL[[x]])/unlist(maxFitL)[x], col=colSet[x], type="l");})
dev.off()

#### Find avg number of drivers to produce a cancer
medFit <- unlist(lapply(bestFitL, median, na.rm=T))
## Recurse through all combinations until their fitness is >= a certain threshold
findNextCancer <- function(start, edgeTab, fit, maxFit) {
  if (fit[start] >= maxFit) {
    return(start)
  }
  offs <- edgeTab$Offspring[which(edgeTab$Parent == start)]
  return(unlist(lapply(offs, findNextCancer, edgeTab, fit, maxFit)))
}

## Not an Endpoint if a parent node's fitness was >= threshold 
isCancerEndPoint <- function(n, edgeTab, fit, maxFit) {
  parents <- edgeTab$Parent[which(edgeTab$Offspring == n)]
  if (length(which(fit[parents] >= maxFit)) > 0) {
    return(0)
  }
  return(1)
}

n0="Node0"
panelSize=15
ngL <- mpsL <- list()
library(beanplot)
library(beeswarm)
pdf(sprintf("plots/network_ndrivers_%d.pdf", panelSize))
for (s in names(bestFitL)) {
  print(s)
  genes <- topL[[s]][1:min(length(topL[[s]]), panelSize)]
  sMat <- read.table(paste("../tcga/clonMat_", s, ".txt", sep=""),
                     sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
  sMat <- sMat[genes,]
  nmut <- apply(sMat, 1, function(x){length(which(x > 0))})
  mpsL[[s]] <- mutPerSamp <- apply(sMat, 2, function(x){length(which(x > 0))})
  ## Load net
  nodeTab <- read.table(sprintf("%s_%d_combFit_nodes.txt",s,panelSize), sep="\t", header=T, stringsAsFactors=F)
  edgeTab <- read.table(sprintf("%s_%d_combFit_edges.txt",s,panelSize), sep="\t", header=T, stringsAsFactors=F)
  fitness <- nodeTab$Fitness
  names(fitness) <- nodeTab$Genes

  fitEnough <- names(fitness)[which(fitness > 1)]
  endIdx <- which(unlist(mclapply(fitEnough, isCancerEndPoint, fit=fitness, edgeTab=edgeTab, maxFit=1, mc.cores=10)) == 1)
  endPoints <- fitEnough[endIdx]
  ##medCombs <- findNextCancer(start="Node0", edgeTab=edgeTab, fit=fitness, maxFit=medFit[s])
  ngL[[s]] <- ngenes <- unlist(lapply(endPoints, function(x){length(unlist(strsplit(x,":")))}))
  wt <- wilcox.test(mutPerSamp, ngenes)
  boxplot(mutPerSamp, ngenes, ylab="Number of drivers", names=c("TCGA", "Landscape"), outline=F,
          ylim=c(0,max(c(mutPerSamp,ngenes))), sub=sprintf("p=%.1e (Wilcoxon test)", wt$p.value))
  beeswarm(c(mutPerSamp, ngenes) ~ rep(1:2, c(length(mutPerSamp), length(ngenes))), pch=pchSet[s], col=colSet[s],
           corral="wrap", add=T)
}

mmps <- unlist(lapply(mpsL, mean))
mng <- unlist(lapply(ngL, mean))
ct <- cor.test(mmps,mng)
mlm <- lm(mng ~ mmps)
plot(mmps, mng, pch=pchSet[names(mpsL)], col=colSet[names(mpsL)], cex=2,
     ylim=c(min(c(unlist(mmps), unlist(mng))), max(c(unlist(mmps), unlist(mng)))),
     xlim=c(min(c(unlist(mmps), unlist(mng))), max(c(unlist(mmps), unlist(mng)))),
     ylab="Mean number of drivers (TCGA)", xlab="Mean number of drivers (Landscapes)",
     sub=sprintf("p=%.1e (R^2=%.2f)", ct$p.value, ct$estimate^2))
abline(mlm, lty=2)

par(mfrow=c(2,1))
beanplot(mpsL, ylab="Number of drivers (TCGA)", outline=F, cut=0.1, las=2, bw=.5, col=rainbow11[c(8,3,4,1)])
## boxplot(mpsL, ylab="Number of drivers (TCGA)", outline=F, ylim=c(0,max(unlist(mpsL))))
## beeswarm(mpsL, pwpch=rep(pchSet[names(mpsL)], unlist(lapply(mpsL, length))),
##          pwcol=rep(colSet[names(mpsL)], unlist(lapply(mpsL, length))), corral="wrap", add=T)

beanplot(ngL, ylab="Number of drivers (Landscapes)", outline=F, cut=0.5, las=2, bw=.5, col=rainbow11[c(8,3,4,1)])
## boxplot(ngL, ylab="Number of drivers (Landscapes)", outline=F, ylim=c(0,max(unlist(ngL))))
## beeswarm(ngL, pwpch=rep(pchSet[names(ngL)], unlist(lapply(ngL, length))),
##          pwcol=rep(colSet[names(ngL)], unlist(lapply(ngL, length))), corral="wrap", add=T)

dev.off()

