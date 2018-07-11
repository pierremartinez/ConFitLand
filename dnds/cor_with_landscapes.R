library(gdata)

load("E:/data/natsuki/selectMutL.RData")
load("E:/data/natsuki/topL.RData")
zapDat <- read.xls("genes_dnds.xlsx")

zapSel <- zapDat$dnds_single
names(zapSel) <- zapDat$Hugo_symbol

rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")
sets=names(topL)
pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- sets
colSet <- rainbow10[1:length(sets)]
names(colSet) <- sets

pdf("plots/landscapes_zapata_cor.pdf")
selMat <- matrix(nr=0, nc=4)
for (s in names(selectMutL)) {
  ov <- intersect(names(selectMutL[[s]]), zapDat$Hugo_symbol)
  ov2 <- intersect(topL[[s]], zapDat$Hugo_symbol)
  selMat <- rbind(selMat, cbind(s, ov2, selectMutL[[s]][ov2], zapSel[ov2]))
  ct = cor.test(selectMutL[[s]][ov2], zapSel[ov2], method="spearman")
  plot(selectMutL[[s]][ov2], zapSel[ov2], main=s, xlab="Selection score", ylab="Zapata dN/dS",
       pch=pchSet[s], col=colSet[s], sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
  abline(lm(zapSel[ov2] ~ selectMutL[[s]][ov2]), lty=2)
  plot(log2(selectMutL[[s]][ov2]), log2(zapSel[ov2]), main=s, xlab="Selection score (log2)", ylab="Zapata dN/dS (log2)",
       pch=pchSet[s], col=colSet[s], sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
  abline(lm(log2(zapSel[ov2]) ~ log2(selectMutL[[s]][ov2])), lty=2)
}
selMat <- data.frame(Set=as.character(selMat[,1]), Gene=as.character(selMat[,2]),
                     Martinez=as.numeric(selMat[,3]), Zapata=as.numeric(selMat[,4]),
                     stringsAsFactors=F)
driv <- unique(selMat[,"Gene"])

ct = cor.test(selMat$Martinez, selMat$Zapata, method="spearman")
plot(selMat$Martinez, selMat$Zapata, main="All sets", xlab="Selection score", ylab="Zapata dN/dS",
     pch=pchSet[selMat$Set], col=colSet[selMat$Set], sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
abline(lm(selMat$Zapata ~ selMat$Martinez), lty=2)
legend("topright", legend=names(pchSet), col=colSet[names(pchSet)], pch=pchSet, bty="n")
plot(log2(selMat$Martinez), log2(selMat$Zapata), main="All sets", xlab="Selection score", ylab="Zapata dN/dS",
     pch=pchSet[selMat$Set], col=colSet[selMat$Set], sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
abline(lm(log2(selMat$Zapata) ~ log2(selMat$Martinez)), lty=2)
legend("topleft", legend=names(pchSet), col=colSet[names(pchSet)], pch=pchSet, bty="n")

meanSel <- unlist(lapply(driv, function(g, selMat) {mean(selMat$Martinez[which(selMat$Gene == g)])}, selMat=selMat))
medSel <- unlist(lapply(driv, function(g, selMat) {median(selMat$Martinez[which(selMat$Gene == g)])}, selMat=selMat))
names(meanSel) <- driv

ct = cor.test(meanSel, zapSel[driv], method="spearman")
plot(meanSel, zapSel[driv], main="Pan cancer drivers", xlab="Mean selection score", ylab="Zapata dN/dS",
     pch=16, col="black", sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
abline(lm(zapSel[driv] ~ meanSel), lty=2, col="red")

plot(log2(meanSel), log2(zapSel[driv]), main="Pan cancer drivers", xlab="Mean selection score (log2)", ylab="Zapata dN/dS (log2)",
     pch=16, col="black", sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
abline(lm(log2(zapSel[driv]) ~ log2(meanSel)), lty=2, col="red")
dev.off()

allg <- alls <-c()
clonL <- list()
for (s in names(topL)) {
  print(s)
  clonMat <- read.table(paste("E:/data/natsuki/clonMat_", s, ".txt", sep=""), sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
  allg <- unique(c(allg, rownames(clonMat)))
  alls <- unique(c(alls, colnames(clonMat)))
  clonL[[s]] <- clonMat
}
allMat <- matrix(0, nr=length(allg), nc=length(alls))
rownames(allMat) <- sort(allg)
colnames(allMat) <- alls
for (s in names(topL)) {
  print(s)
  clonMat <- clonL[[s]]
  allMat[rownames(clonMat), colnames(clonMat)][clonMat>0] <- clonMat[clonMat>0]
}
save(allMat, file="allMat.RData")

##rm(list=c("allMat", "clonMat", "clonL"))

getExpectMut <- function(g, totMut, len) {
  (len[g] / sum(len, na.rm=T)) * totMut
}

## gen.ann: Gene annotations (transcript length in particular)
load("E:/data/genes_annotation.RData")
glen <- gen.ann$median_trancript_length
names(glen) <- gen.ann$gene_symbol
## hack for wrong symbols
glen[c("FAM123B", "HNRPDL", "MLL2", "MLL3", "EIF2C3", "MLL", "CSDA", "C15orf55")] <-
  glen[c("AMER1", "HNRNPDL", "KMT2D", "KMT2C", "AGO3", "KMT2A", "YBX3", "NUTM1")]
## Data on select CN events
load("E:/data/TCGA/cn_drivers/cnSelectL.RData")

nmut <- apply(allMat, 1, function(x){length(which(x > 0))})
expectMut <- unlist(lapply(rownames(allMat)[grep("_", rownames(allMat), invert=T)],
                           getExpectMut, totMut=sum(nmut), len=glen[rownames(allMat)]))
selectAll <- nmut[grep("_", rownames(allMat), invert=T)] / expectMut
##selectMut <- c(selectMutL[[s]], cnSelectL[[s]][rownames(allMat)[grep("_", rownames(allMat))]])

pdf("plots/landscapes_zapata_pancancer_cor.pdf")
ov3 <- intersect(zapDat$Hugo_symbol, names(selectAll))
ct = cor.test(selectAll[ov3], zapSel[ov3], method="spearman")
plot(selectAll[ov3], zapSel[ov3], main=s, xlab="Selection score", ylab="Zapata dN/dS",
     pch=16, col="black", sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
abline(lm(zapSel[ov3] ~ selectAll[ov3]), lty=2)
plot(log2(selectAll[ov3]), log2(zapSel[ov3]), main=s, xlab="Selection score (log2)", ylab="Zapata dN/dS (log2)",
     pch=16, col="black", sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
abline(lm(log2(zapSel[ov3]) ~ log2(selectAll[ov3])), lty=2)

ct = cor.test(selectAll[driv], zapSel[driv], method="spearman")
plot(selectAll[driv], zapSel[driv], main=s, xlab="Selection score", ylab="Zapata dN/dS",
     pch=16, col="black", sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
abline(lm(zapSel[driv] ~ selectAll[driv]), lty=2)
plot(log2(selectAll[driv]), log2(zapSel[driv]), main="Pan cancer drivers", xlab="Selection score (log2)", ylab="Zapata dN/dS (log2)",
     pch=16, col="black", sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
abline(lm(log2(zapSel[driv]) ~ log2(selectAll[driv])), lty=2, col="red")
dev.off()


