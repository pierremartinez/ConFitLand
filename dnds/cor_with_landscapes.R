library(gdata)

load("../tcga/selectMutL.RData")
load("../tcga/topL.RData")
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

selMat <- matrix(nr=0, nc=4)
for (s in names(selectMutL)) {
  ov <- intersect(names(selectMutL[[s]]), zapDat$Hugo_symbol)
  ov2 <- intersect(topL[[s]], zapDat$Hugo_symbol)
  selMat <- rbind(selMat, cbind(s, ov2, selectMutL[[s]][ov2], zapSel[ov2]))
}
selMat <- data.frame(Set=as.character(selMat[,1]), Gene=as.character(selMat[,2]),
                     Martinez=as.numeric(selMat[,3]), Zapata=as.numeric(selMat[,4]),
                     stringsAsFactors=F)
driv <- unique(selMat[,"Gene"])

pdf("plots/landscapes_zapata_cor.pdf")
ct = cor.test(selMat$Martinez, selMat$Zapata, method="spearman")
plot(log2(selMat$Martinez), log2(selMat$Zapata), main="All sets", xlab="Selection score", ylab="Zapata dN/dS",
     pch=pchSet[selMat$Set], col=colSet[selMat$Set], sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
abline(lm(log2(selMat$Zapata) ~ log2(selMat$Martinez)), lty=2)
legend("topleft", legend=names(pchSet), col=colSet[names(pchSet)], pch=pchSet, bty="n")

meanSel <- unlist(lapply(driv, function(g, selMat) {mean(selMat$Martinez[which(selMat$Gene == g)])}, selMat=selMat))
medSel <- unlist(lapply(driv, function(g, selMat) {median(selMat$Martinez[which(selMat$Gene == g)])}, selMat=selMat))
names(meanSel) <- driv

ct = cor.test(meanSel, zapSel[driv], method="spearman")
plot(log2(meanSel), log2(zapSel[driv]), main="Pan cancer drivers", xlab="Mean selection score (log2)", ylab="Zapata dN/dS (log2)",
     pch=16, col="black", sub=sprintf("rho=%.2f (p=%.1e)", ct$estimate, ct$p.value))
abline(lm(log2(zapSel[driv]) ~ log2(meanSel)), lty=2, col="red")
dev.off()



