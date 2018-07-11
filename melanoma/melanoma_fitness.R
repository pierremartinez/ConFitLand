## Melanoma & nevi fitness
## Be careful to adapt the getFitStr function to the appropriate one if using one of the separate models.

library(gdata)
## Load mutational data
mutDat <- read.xls("nejmoa1502583_appendix_4.xlsx", skip=1, stringsAsFactors=F)
mutDat$Histology.1 <- gsub("Desmplastic", "Desmoplastic", mutDat$Histology.1)
## Reformatting histology
allHist <- c("Nevus",
             "IntermediateMalignant",
             "Melanoma in situ",
             "Melanoma in situ - parts 1 and 2",
             "IntermediateBenign",
             "Melanoma in situ - part 1",
             "Blue Nevus",
             "Melanoma(T1)",
             "Melanoma(T2)",
             "Melanoma in situ - parts 3 and 4",
             "Desmoplastic Melanoma",
             "Melanoma in situ - part 2",
             "Blue Nevus Like Melanoma",
             "Melanoma(T3)",
             "Melanoma(T4)",
             "Metastasis")
equivHist <- c("Nevus",
             "IntermediateMalignant",
             "MIS",
             "MIS",
             "IntermediateBenign",
             "MIS",
             "BN",
             "Melanoma(T1)",
             "Melanoma(T2)",
             "MIS",
             "DesmoplasticMelanoma",
             "MIS",
             "BNLM",
             "Melanoma(T3)",
             "Melanoma(T4)",
             "Metastasis")
names(equivHist) <- allHist

## Which annotation is benign / malign
benHist <- c("Nevus",
             "BN",
             "Blue Nevus",
             "IntermediateBenign")
malHist <- c("IntermediateMalignant",
             "MIS",
             "BNLM",
             "DesmoplasticMelanoma",
             "Melanoma(T1)",
             "Melanoma(T2)",
             "Melanoma(T3)",
             "Melanoma(T4)")

## Score sample malignancy
malignState <- c(0,0,0,1, 2,3,4,4,4,5,6,7)
names(malignState) <- c(benHist, malHist)

## Load CNA data
cnvDat <- read.xls("nejmoa1502583_appendix_5.xlsx", stringsAsFactors=F)
cnvPat <- unlist(lapply(cnvDat$ID, function(x){unlist(strsplit(as.character(x),split="_"))[1]}))
cnvType <- unlist(lapply(cnvDat$ID, function(x){unlist(strsplit(as.character(x),split="_"))[2]}))
cnvN <- unlist(lapply(cnvDat$ID, function(x){unlist(strsplit(as.character(x),split="_"))[3]}))

normCN <- cnvDat$seg.mean[which(cnvType=="Normal")]
normCN <- normCN[which(abs(normCN) < 2)]
mCN <- mean(normCN, na.rm=T)
sdCN <- sd(normCN, na.rm=T)

#### Get IntOGen drivers
sets <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "SKCM")
intomut <- read.table("../misc/Mutational_drivers_per_tumor_type.tsv", sep="\t", header=T, stringsAsFactors=F, check.names=F)
intocna <- read.table("../misc/CNA_drivers_per_tumor_type.tsv", sep="\t", header=T, stringsAsFactors=F, check.names=F)
equivSet <- sets
names(equivSet) <- sets
equivSet["COAD"] <- "COREAD"
equivSet["SKCM"] <- "CM"
equivSet["KIRC"] <- "RCCC"
drivL <- list()
for (s in sets) {
  tmp <- paste(intocna$geneHGNCsymbol, intocna$GISTIC_CNA, sep="_")[which(intocna$Tumor_type_GISTIC == equivSet[s])]
  tmp <- gsub("_A$", "_gain", gsub("_D$", "_loss", tmp))
  drivL[[s]] <- c(intomut$geneHGNCsymbol[which(intomut$Tumor_type == equivSet[s])], tmp)
}

cmcn <- drivL[["SKCM"]][grep("_", drivL[["SKCM"]])]
cneg <- unique(unlist(lapply(cmcn, function(x){unlist(strsplit(x, "_"))[1]})))
## Get driver coordinates 
library(biomaRt)
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)
ensres <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                filter="hgnc_symbol", values=cneg, mart=ensembl)

pats <- unique(cnvPat)
patType <- paste(cnvPat, cnvType, sep="_")
validIds <- unique(patType[grep("Normal", patType, invert=T)])
repPerId <- unlist(lapply(validIds, function(id) {
  idid <- which(paste(cnvPat, cnvType, sep="_") == id)
  return(length(unique(cnvN[idid])))
}))
names(repPerId) <- validIds

## Investigate CNA presence in all samples
cnMat <- matrix(0, nr=length(cmcn), nc=length(validIds))
colnames(cnMat) <- validIds
rownames(cnMat) <- cmcn
for (cn in cmcn) {
  l <- unlist(strsplit(cn, split="_"))
  st <- as.numeric(ensres[which(ensres$hgnc_symbol == l[1]),"start_position"])
  e <- as.numeric(ensres[which(ensres$hgnc_symbol == l[1]),"end_position"])
  chr <- ensres[which(ensres$hgnc_symbol == l[1]),"chromosome_name"]
  if (chr == "X")
    chr=23
  if (chr == "Y")
    chr=24
  chr=as.numeric(chr)
  idx1 <- which(cnvDat$chrom == chr & patType %in% validIds)
  subcnvDat <- cnvDat[idx1,]
  idx2 <- which(subcnvDat$loc.start < e & subcnvDat$loc.end > st)
  subcnvDat <- subcnvDat[idx2,]
  subIds <- paste(cnvPat, cnvType, sep="_")[idx1[idx2]]
  
  lt <- ifelse(l[2] == "gain", FALSE, TRUE)
  pvals <- unlist(lapply(subcnvDat$seg.mean, pnorm, m=mCN, sd=sdCN, lower.tail=lt))
  
  for (id in subIds[which(pvals < 0.05)])
    cnMat[cn,id] <- cnMat[cn,id] + (1 / repPerId[id])
}
cnMat[cnMat<1] <- 0
write.table(cnMat, file="Bastian_SKCM_cnDriv_table.txt", sep="\t", quote=F, col.names=NA)

## Minimum value for mutations to be clonal
clonMAF <- .8
mutIds <- unique(c(paste(mutDat$Sample, equivHist[as.character(mutDat$Histology)], sep="_"),
                   paste(mutDat$Sample, equivHist[as.character(mutDat$Histology.1)], sep="_"),
                   paste(mutDat$Sample, equivHist[as.character(mutDat$Histology.2)], sep="_")))
genes <- unique(mutDat$Hugo.Symbol)
## Build clonal alteration matrices
mutMat <- matrix(0, nr=length(genes), nc=length(mutIds))
colnames(mutMat) <- mutIds
rownames(mutMat) <- genes
for (i in 1:nrow(mutDat)) {
  h <- equivHist[as.character(mutDat$Histology[i])]
  p <- mutDat$Sample[i]
  id <- paste(p, h, sep="_")
  if (mutDat$NormalizedMAF[i] >= clonMAF)
    mutMat[mutDat$Hugo.Symbol[i],id] <- 1
  if (mutDat$Histology.1[i] != "") {
    h1 <- equivHist[as.character(mutDat$Histology.1[i])]
    id1 <- paste(p, h1, sep="_")
    if (mutDat$NormalizedMAF.1[i] >= clonMAF)
      mutMat[mutDat$Hugo.Symbol[i],id1] <- 1
  }
  if (mutDat$Histology.2[i] != "") {
    h2 <- equivHist[as.character(mutDat$Histology.2[i])]
    id2 <- paste(p, h2, sep="_")
    if (mutDat$NormalizedMAF.2[i] >= clonMAF)
      mutMat[mutDat$Hugo.Symbol[i],id2] <- 1
  }
}
write.table(mutMat, file="Bastian_SKCM_mutation_table.txt", sep="\t", quote=F, col.names=NA)

## Indices to know which sample/patient is which
clonIds <- intersect(validIds, mutIds)
clonHist <- unlist(lapply(clonIds, function(x){unlist(strsplit(x, "_"))[2]}))
clonPats <- unlist(lapply(clonIds, function(x){unlist(strsplit(x, "_"))[1]}))
clonMalign <- unlist(lapply(clonHist, function(h){h %in% malHist}))
clonBenign <- unlist(lapply(clonHist, function(h){h %in% benHist}))
names(clonMalign) <- names(clonBenign) <- clonIds
validIdx <- which(clonBenign | clonMalign)
validClon <- clonIds[validIdx]
validPats <- clonPats[validIdx]
validHist <- clonHist[validIdx]
validMalign <- clonMalign[validIdx]
validBenign <- clonBenign[validIdx]

clonMat <- rbind(mutMat[,validClon], cnMat[,validClon])
write.table(clonMat, file="Bastian_SKCM_clonMat.txt", sep="\t", quote=F, col.names=NA)

#############################
####   Compute Fitness   ####
#############################

## Update if needed
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
load("../misc/cnSelectL.RData")

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

s="SKCM"
## Update if using separate model
## occnf = bestCombs[s,1]
## socnf = bestCombs[s,2]
## epinf = bestCombs[s,3]
intnf = bestCombs[s,1]
epinf = bestCombs[s,2]
maxland = 50
medFit <- unlist(lapply(bestFitL, median, na.rm=T))
w <- rep(1/unlist(lapply(bestFitL, length)), unlist(lapply(bestFitL, length)))
medFit2 <- weighted.mean(unlist(bestFitL), w=w, na.rm=T)
names(medFit) <- names(bestFitL)
fitMat <- matrix(NA, nr=length(unique(validPats)), nc=2)
rownames(fitMat) <- unique(validPats)
colnames(fitMat) <- c("benFit", "malFit")

library(parallel)
library(beeswarm)
rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")
genes <- intersect(topL[[s]][1:min(length(topL[[s]]), maxland)], rownames(clonMat))
sMat <- clonMat[genes,]
skcmMat <- read.table(paste("../tcga/clonMat_", s, ".txt", sep=""),
                      sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
nmut <- apply(skcmMat, 1, sum)
selectMut <- c(selectMutL[[s]], cnSelectL[[s]][rownames(skcmMat)[grep("_", rownames(skcmMat))]])
sMat <- clonMat[genes,]
mutPerSamp <- apply(sMat, 2, sum)
socRat <- (1 / socRatL[[s]][genes]) ^ socPowL[[s]][genes]
names(socRat) <- genes
fMat <- read.table(paste("../tcga/epiFactMat_", s, ".txt", sep=""), sep="\t", header=T, row.names=1, check.names=F)
fMat <- fMat[genes,]
geneComb <- apply(sMat, 2, function(x){paste(rownames(sMat)[which(x==1)], collapse=":")})
## Update if using separate model
## sFit <- unlist(mclapply(geneComb, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, of=occnf, sf=socnf, ef=epinf, mc.cores=8))
sFit <- unlist(mclapply(geneComb, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, sf=intnf, ef=epinf, mc.cores=8))
sFitN <- sFit / median(medFit[s], na.rm=T)

for (p in unique(validPats)) {
  pben <- validClon[intersect(grep(paste("^",p,"_",sep=""), validClon), which(validBenign))]
  if (length(pben) > 1) {
    print(paste("Error", p))
  }
  pmal <- validClon[intersect(grep(paste("^",p,"_",sep=""), validClon), which(validMalign))]
  if (length(pmal) > 1) {
    ## Select least advanced sample
    m <- malignState[unlist(lapply(pmal, function(x){unlist(strsplit(x,"_"))[2]}))]
    pmal=pmal[which(m == min(m))[1]]
  }
  if (length(pben) > 0) {
    fitMat[p,"benFit"] <- sFitN[pben]
  }
  if (length(pmal) > 0) {
    fitMat[p,"malFit"] <- sFitN[pmal]
  }
}

pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- sets
colSet <- rainbow10[1:length(sets)]
names(colSet) <- sets

## Unpaired version
## bothIdx <- which(!is.na(fitMat[,"benFit"]) & !is.na(fitMat[,"malFit"]))
## ptt <- t.test(fitMat[bothIdx,"benFit"], fitMat[bothIdx,"malFit"], paired=T)
## wt <- wilcox.test(fitMat[bothIdx,"benFit"], fitMat[bothIdx,"malFit"])
## wt2 <- wilcox.test(fitMat[,"benFit"], fitMat[,"malFit"])

## pdf(paste("plots/bastian_melanoma_",model,"_fitness.pdf",sep=""))
## boxplot(fitMat[bothIdx,"malFit"], fitMat[bothIdx,"benFit"], ylab="Normalised fitness", names=c("Malign", "Benign"),
##         ylim=c(0, max(fitMat[bothIdx,])), outline=F, sub=sprintf("p=%.3f (paired t-test)", ptt$p.value))
## beeswarm(c(fitMat[bothIdx,"malFit"], fitMat[bothIdx,"benFit"]) ~ rep(1:2, rep(length(bothIdx), 2)),
##          pch=pchSet[s], pwcol=rep(rainbow10[c(1,7)], rep(length(bothIdx), 2)),
##          add=T, cex=2, corral="wrap")

## boxplot(fitMat[,"malFit"], fitMat[,"benFit"], ylab="Normalised fitness", names=c("Malign", "Benign"),
##         ylim=c(0, max(fitMat, na.rm=T)), outline=F, sub=sprintf("p=%.3f (Wilcoxon test)", wt2$p.value))
## beeswarm(c(fitMat[,"malFit"], fitMat[,"benFit"]) ~ rep(1:2, rep(nrow(fitMat),2)),
##          pch=pchSet[s], pwcol=rep(rainbow10[c(1,7)], rep(nrow(fitMat), 2)), add=T, cex=2, corral="wrap")
## dev.off()

precMat <- matrix(NA, nr=length(unique(validPats)), nc=2)
rownames(precMat) <- unique(validPats)
colnames(precMat) <- c("benFit", "malFit")
for (p in unique(validPats)) {
  pall <- validClon[grep(paste("^",p,"_",sep=""), validClon)]
  
  if (length(pall) <= 1) {
    if (pall %in% validClon[validBenign]) {
      
    }
  } else {
    m <- malignState[unlist(lapply(pall, function(x){unlist(strsplit(x,"_"))[2]}))]
    pmal=pall[which(m == max(m))[1]]
    pben=pall[which(m == min(m))[1]]
    precMat[p,"benFit"] <- sFitN[pben]
    precMat[p,"malFit"] <- sFitN[pmal]
  }
}
ptt <- t.test(precMat[,"benFit"], precMat[,"malFit"], paired=T)
wt <- wilcox.test(precMat[,"benFit"], precMat[,"malFit"], paired=T)

pdf(paste("plots/bastian_precursor_fitness_",model,".pdf",sep=""))
boxplot(precMat[,"malFit"], precMat[,"benFit"], ylab="Normalised fitness", names=c("Malign", "Precursor"),
        ylim=c(0, max(precMat,na.rm=T)), outline=F, sub=sprintf("p=%.3f (paired Wilcoxon test)", wt$p.value))
beeswarm(c(precMat[,"malFit"], precMat[,"benFit"]) ~ rep(1:2, rep(nrow(precMat), 2)),
         pch=pchSet[s], pwcol=rep(rainbow10[c(1,7)], rep(nrow(precMat), 2)),
         add=T, cex=2, corral="wrap")
dev.off()

getExpectMut <- function(g, totMut, len) {
  (len[g] / sum(len, na.rm=T)) * totMut
}

## gen.ann: Gene annotations (transcript length in particular)
load("../misc/genes_annotation.RData")
glen <- gen.ann$median_transcript_length
names(glen) <- gen.ann$gene_symbol
## hack for wrong symbols
glen[c("FAM123B", "HNRPDL", "MLL2", "MLL3", "EIF2C3", "MLL", "CSDA", "C15orf55")] <-
  glen[c("AMER1", "HNRNPDL", "KMT2D", "KMT2C", "AGO3", "KMT2A", "YBX3", "NUTM1")]

library(parallel)

## Fitness in other sets
nevMat <- clonMat[,validBenign]
s="SKCM"
nevFitN <- sFitN[validBenign]
osets <- setdiff(names(bestFitNL), s)
ofitMat <- matrix(nr=ncol(nevMat), nc=length(osets))
rownames(ofitMat) <- colnames(nevMat)
colnames(ofitMat) <- c(osets)
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
  expectMut <- unlist(lapply(rownames(sMat)[grep("_", rownames(sMat), invert=T)],
                             getExpectMut, totMut=sum(nmut), len=glen[rownames(sMat)]))
  selectMut <- nmut[grep("_", rownames(sMat), invert=T)] / expectMut
  selectMut <- c(selectMut, cnSelectL[[s]][rownames(sMat)[grep("_", rownames(sMat))]])
  genes <- intersect(genes, rownames(nevMat))
  socRat <- (1 / socRatL[[s]][genes]) ^ socPowL[[s]][genes]
  names(socRat) <- genes

  ssMat <- nevMat[genes,]
  mutPerSamp <- apply(ssMat, 2, sum)
  geneComb <- apply(ssMat, 2, function(x){paste(rownames(ssMat)[which(x==1)], collapse=":")})
  oFit <- unlist(mclapply(geneComb, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, sf=ssf, ef=sef, mc.cores=8))
  oFitN <- oFit / median(medFit[s], na.rm=T)
  ofitMat[,s] <- oFitN
}

medOther <- apply(ofitMat, 1, median)
meanOther <- apply(ofitMat, 1, mean)

s="SKCM"
pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- names(socPowL)
colSet <- rainbow10[1:length(socPowL)]
names(colSet) <- names(socPowL)
pdf(paste("plots/nevi_other_fitness_",model,".pdf", sep=""))
par(mar=c(6,4,1,1))
boxplot(t(cbind(nevFitN, ofitMat)), ylim=c(0, max(c(nevFitN, ofitMat))), las=2, ylab="Normalised Fitness", outline=F)
beeswarm(c(nevFitN, ofitMat) ~ rep(1:ncol(nevMat), length(osets)+1), add=T, corral="wrap",
         pwpch=rep(pchSet[c(s,osets)], rep(nrow(ofitMat), length(osets)+1)),
         pwcol=rep(colSet[c(s,osets)], rep(nrow(ofitMat), length(osets)+1)))
legend("topright", legend=names(bestFitNL), col=colSet[names(bestFitNL)], pch=pchSet[names(bestFitNL)])

ptt <- t.test(nevFitN, medOther, paired=T)
boxplot(cbind(nevFitN, medOther), outline=F, ylim=c(0, max(c(nevFitN, medOther))), names=c(s, "Other sets (median)"),
        ylab="Normalised Fitness", sub=sprintf("p=%.1e (paired t-test)", ptt$p.value))
beeswarm(c(nevFitN, medOther) ~ rep(1:2, rep(length(nevFitN), 2)), pch=pchSet[s], col=colSet[s], add=T, cex=2)

ptt <- t.test(nevFitN, meanOther, paired=T)
boxplot(cbind(nevFitN, meanOther), outline=F, ylim=c(0, max(c(nevFitN, meanOther))), names=c(s, "Other sets (mean)"),
        ylab="Normalised Fitness", sub=sprintf("p=%.1e (paired t-test)", ptt$p.value))
beeswarm(c(nevFitN, meanOther) ~ rep(1:2, rep(length(nevFitN), 2)), pch=pchSet[s], col=colSet[s], add=T, cex=2)
dev.off()
