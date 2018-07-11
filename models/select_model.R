allmodels <- c("separate_mean",
               "separate_prod",
               "combined_mean",
               "combined_prod")

## Most clonally mutated genes
load("../tcga/topL.RData")
maxland=25
ctL <- ctpL <- sdL <- spearL <- spearpL <- list()
for (m in allmodels) {
  print(m)
  load(paste("bestCombs_",m,".RData", sep=""))
  load(paste("bestFitNL_",m,".RData", sep=""))

  ctL[[m]] <- c()
  for (s in names(bestFitNL)) {
    print(s)
    sFitN <- bestFitNL[[s]]
    sdL[[m]][s] <- sd(sFitN)
    genes <- topL[[s]][1:min(length(topL[[s]]), maxland)]
    clonMat <- read.table(paste("../clonMat_", s, ".txt", sep=""), sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
    drivPerSamp <- apply(clonMat[genes,], 2, sum)
    ctL[[m]][s] <- cor.test(drivPerSamp, sFitN)$estimate ^ 2
    ctpL[[m]][s] <- cor.test(drivPerSamp, sFitN)$p.value
    spearL[[m]][s] <- cor.test(drivPerSamp, sFitN, method="spearman")$estimate
    spearpL[[m]][s] <- cor.test(drivPerSamp, sFitN, method="spearman")$p.value
  }
}

library(beeswarm)
rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")
pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- names(bestFitNL)
colSet <- rainbow10[1:length(bestFitNL)]
names(colSet) <- names(bestFitNL)

pchM <- c(7:14)[1:length(ctL)]
names(pchM) <- names(ctL)
colM <- rainbow11[1:length(ctL)]
names(pchM) <- names(ctL)

pdf("plots/model_selection.pdf")
par(mar=c(10,4,.5,.5))
boxplot(ctL, las=2, ylab="Correlation between fitness and number of drivers (R^2)", outline=F, ylim=c(0,1))
beeswarm(unlist(ctL) ~ rep(1:length(ctL), unlist(lapply(ctL, length))), add=T,
         pwpch=pchSet[unlist(lapply(ctL, names))], pwcol=colSet[unlist(lapply(ctL, names))])
boxplot(spearL, las=2, ylab="Correlation between fitness and number of drivers (Spearman's rho)", outline=F, ylim=c(0,1))
beeswarm(unlist(spearL) ~ rep(1:length(spearL), unlist(lapply(spearL, length))), add=T,
         pwpch=pchSet[unlist(lapply(spearL, names))], pwcol=colSet[unlist(lapply(spearL, names))])
boxplot(sdL, las=2, ylab="Standard deviation", outline=F, ylim=c(min(unlist(sdL)), max(unlist(sdL))))
beeswarm(unlist(sdL) ~ rep(1:length(sdL), unlist(lapply(sdL, length))), add=T,
         pwpch=pchSet[unlist(lapply(sdL, names))], pwcol=colSet[unlist(lapply(sdL, names))])
medCt <- unlist(lapply(ctL, median))
medSpear <- unlist(lapply(spearL, median))
medSd <- unlist(lapply(sdL, median))
plot(medCt, medSd, las=2, ylab="Median standard deviation", xlab="Median R^2", pch=pchM, cex=2, col=colM)
legend("bottomleft", pch=pchM, legend=names(pchM), col=colM)
plot(medSpear, medSd, las=2, ylab="Median standard deviation", xlab="Median rho", pch=pchM, cex=2, col=colM)
legend("bottomleft", pch=pchM, legend=names(pchM), col=colM)
dev.off()

