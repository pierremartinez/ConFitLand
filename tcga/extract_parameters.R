## Script 1
## Extracting the parameters that will be used to calculate fitness
## Careful to mclapply calls, change mc.cores parameter if not enough processors are available.

library(pwr)
library(plotrix)
library(parallel)
library(beeswarm)
library(beanplot)
alpha = 0.05 ## prob type I error
minPercDriv <- 0.005 ## Mininmum percentage of cohort harbouring each driver

## Plot colouring
rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")

#### Get IntOGen drivers
sets <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "SKCM")
intomut <- read.table("../misc/Mutational_drivers_per_tumor_type.tsv",
                   sep="\t", header=T, stringsAsFactors=F, check.names=F)
intocna <- read.table("../misc/CNA_drivers_per_tumor_type.tsv",
                   sep="\t", header=T, stringsAsFactors=F, check.names=F)
equivSet <- c("BLCA", "BRCA", "COREAD", "GBM", "HNSC", "RCCC", "LUAD", "LUSC", "CM")
names(equivSet) <- sets
drivL <- list()
for (s in sets) {
  tmp <- paste(intocna$geneHGNCsymbol, intocna$GISTIC_CNA, sep="_")[which(intocna$Tumor_type_GISTIC == equivSet[s])]
  tmp <- gsub("_A$", "_gain", gsub("_D$", "_loss", tmp))
  drivL[[s]] <- c(intomut$geneHGNCsymbol[which(intomut$Tumor_type == equivSet[s])], tmp)
}

## For each event e, get number of partners distribution
getNPartners <- function(e, tab) {
  idxe <- which(tab[e,] != 0)
  others <- setdiff(rownames(tab), e)
  return(apply(tab[others, idxe], 2, function(x){length(which(x != 0))}))
}

## Estimate power of statistical test
comparePower <- function(i, l) {
  v <- unlist(lapply(setdiff(1:length(l), i), function(j){l[j]}))
  tt <- t.test(l[[i]], v)
  ptt <- power.t.test(n=length(l[[i]]), delta=mean(l[[i]]) - mean(v), sd=sd(v), type="two.sample")
  return(ptt$power)
}

## ratio of mean observed n partners compared to mean expected n partners
partnerRatio <- function(i, l) {
  v <- unlist(lapply(setdiff(1:length(l), i), function(j){l[j]}))
  m <- mean(l[[i]])
  if (m == 0) {
    print(paste(candidates[i], "has no partners "))
    m = 1 / (length(l[[i]]) + 1)
  }
  return(m / mean(v))
}

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
## Data on select CN events
load("cnSelectL.RData")

pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- sets
colSet <- rainbow10[1:length(sets)]
names(colSet) <- sets

## Wilcox test: distribution of number of partners of a given driver vs all other drivers
comparePartners <- function(i, l) {
  v <- unlist(lapply(setdiff(1:length(l), i), function(j){l[j]}))
  return(wilcox.test(l[[i]], v)$p.value)
}

#### Get sociability index
socRatL <- list() ## list of ratio of mean observed n partners compared to mean expected n partners
selectMutL <- list() ## list of ratio of mean observed n partners compared to mean expected n partners
socPowL <- list() ## list of power values
socPL <- list() ## list of p-values
candL <- list() ## list of candidates
maxcand <- 100 ## Maximum number of drivers
maxdisp <- 20 ## Maximum number of drivers displayed in plots
library(beanplot)
pdf("plots/sociability_plot.pdf")
for (i in 1:length(sets)) {
  s=sets[i]
  print(s)
  driv <- drivL[[s]]
  clonMat <- read.table(paste("clonMat_", s, ".txt", sep=""), sep="\t", row.names=1, header=T, stringsAsFactors=F, check.names=F)
  nmut <- apply(clonMat, 1, function(x){length(which(x > 0))})
  expectMut <- unlist(lapply(rownames(clonMat)[grep("_", rownames(clonMat), invert=T)],
                             getExpectMut, totMut=sum(nmut), len=glen[rownames(clonMat)]))
  selectMutL[[s]] <- nmut[grep("_", rownames(clonMat), invert=T)] / expectMut
  selectMut <- c(selectMutL[[s]], cnSelectL[[s]][rownames(clonMat)[grep("_", rownames(clonMat))]])
  mutPerSamp <- apply(clonMat, 2, function(x){length(which(x > 0))})
  candidates <- driv[which(nmut[driv] > max(1, ncol(clonMat) * minPercDriv))]
  drivMat <- clonMat[candidates,]
  candidates <- names(sort(nmut[candidates], decreasing=T))[1:min(length(candidates), maxcand)]
  npartners <- lapply(candidates, getNPartners, tab=drivMat)
  candL[[s]] <-  candidates
  names(npartners) <- candidates

  ## Nb partners & power calculation
  socRatL[[s]] <- unlist(lapply(1:length(npartners), partnerRatio, l=npartners))
  socPowL[[s]] <- unlist(lapply(1:length(npartners), comparePower, l=npartners))
  names(socRatL[[s]]) <- names(socPowL[[s]]) <- names(npartners)

  bestcand <- candidates[1:min(length(candidates), maxdisp)]
  meannp <- unlist(lapply(npartners[bestcand], mean))
  nobs <- unlist(lapply(npartners[bestcand], length))
  names(meannp) <- names(npartners[bestcand])
  meanall <- mean(unlist(npartners[bestcand]))
  npo <- order(meannp)
  wilp <- unlist(lapply(1:length(npartners), comparePartners, l=npartners))
  names(wilp) <- names(npartners)
 
  ## label display tweak
  dispnames=paste(names(npartners[bestcand]), ": ", nobs, sep="")
  
  ## Colour vector
  colv <- rep(colors()[394], length(npartners[bestcand]))
  colv[wilp[npo] < 0.05 & meannp[npo] < meanall] <- "blue"
  colv[wilp[npo] < 0.05 & meannp[npo] > meanall] <- "red"
  
  ## Plotting
  par(mar=c(10,4,1,1))
  beanplot(npartners[bestcand], at=order(npo), what=c(0,1,1,0), ylab="Number of partners", las=2, overallline="median",
           col=c("grey", "black", "white", "darkgreen"), names=dispnames, main=s, bw=.5,
           cutmin=0, cutmax=max(unlist(npartners[bestcand])))
  abline(h=mean(unlist(npartners[bestcand])), lty=3)
  points(sort(meannp), pch=18, col=colv, cex=1.5)
}

allSoc <- unlist(lapply(names(candL), function(s){(1 / socRatL[[s]][candL[[s]]]) ^ socPowL[[s]][candL[[s]]]}))
allSel <- unlist(lapply(names(candL), function(s){selectMutL[[s]][candL[[s]]]}))
allSet <- unlist(lapply(names(candL), function(s){rep(s, length(candL[[s]]))}))
ct <- cor.test(allSoc, allSel)
ct2 <- cor.test(allSoc, log(allSel))
plot(allSoc, allSel, pch=pchSet[allSet], col=colSet[allSet], xlab="Sociability ratio", ylab="Selection score",
     sub=sprintf("R^2=%.2f (p=%1.e)", ct$estimate ^ 2, ct$p.value))
plot(allSoc, log(allSel), pch=pchSet[allSet], col=colSet[allSet], xlab="Sociability ratio", ylab="Selection score (log scale)",
     sub=sprintf("R^2=%.2f (p=%1.e)", ct2$estimate ^ 2, ct2$p.value))
dev.off()

save(candL, file="candL.RData")
save(socRatL, file="socRatL.RData")
save(socPowL, file="socPowL.RData")
save(selectMutL, file="selectMutL.RData")

hyperInSamp <- function(i, g, mutPerSamp, nmut) {
  ## probability that at least one mutation in sample i is gene g (>0, hence the upper tail)
  phyper(m=nmut[g], n=(sum(mutPerSamp) - nmut[g]), k=mutPerSamp[i], q=0, lower.tail=F)
}

drawSuccess <- function(p) {
  sample(x=c(0,1), size=1, prob=c(1-p, p))
}

randDistrib <- function(probv, nsim=1000, mc.cores=1) {
  unlist(mclapply(1:nsim, function(i, probv){sum(unlist(lapply(probv, drawSuccess)))}, probv=probv, mc.cores=mc.cores))
}

nsim=10000

##load("cnSelectL.RData")

getExpectMut <- function(g, totMut, len) {
  (len[g] / sum(len, na.rm=T)) * totMut
}

library(parallel)
freqL <- list() ## List of frequencies of each alteration per set
obsL <- list() ## List of number of occurences of each combination per set
unionL <- list() ## List of samples impacted by any of the two alteration for each combination per set
ratL <- list() ## List of ratios observed / expected per combination per set
drawpL <- list() ## List of p-values estimated from percentage of distrib < or > to observed value
esL <- list() ## List of effect size per combination
rdrawL <- list() ## List of random draws per combination
combL <- list() ## List of gene combinations
topL <- list() ## List of top 3 genes per set
pow2pL <- list() ## List of power achieved
maxtop <- 50
for (i in 1:length(sets)) {
  s=sets[i]
  print(s)
  clonMat <- read.table(paste("clonMat_", s, ".txt", sep=""), sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)

  ## Get most selected alterations
  nmut <- apply(clonMat, 1, sum)
  mutPerSamp <- apply(clonMat, 2, sum)
  candidates <- candL[[s]]
  top <- names(sort(nmut[candidates], decreasing=T))
  ## Filter out CDKN2B loss when CDKN2A and CDKN2B loss are drivers (same chr location).
  tmpa <- grep("CDKN2A_loss", top)
  tmpb <- grep("CDKN2B_loss", top)
  if (length(tmpa) == 1 & length(tmpb) == 1)
    top <- top[-tmpb]
  top <- top[which(!is.na(top))]
  top <- top[1:min(maxtop,length(top)-1)]
  topL[[s]] <- top
  freqL[[s]] <- nmut[top]/ncol(clonMat)

  ## Get pairwise combinations
  combs <- combn(length(top),2)
  rdrawL[[s]] <- list()
  ratL[[s]] <- obsL[[s]] <- unionL[[s]] <- drawpL[[s]] <- esL[[s]] <- combL[[s]] <- c()
  pow2pL[[s]] <- list()
  for (g in top)
    pow2pL[[s]][[g]] <- c()
  for (j in 1:ncol(combs)) {
    ## Get expectations of co-occurrences
    g1=top[combs[1,j]]
    g2=top[combs[2,j]]
    combL[[s]] <- c(combL[[s]], paste(g1,"-",g2))

    n1 <- length(which(clonMat[g1,] > 0))
    n2 <- length(which(clonMat[g2,] > 0))
    print.noquote(paste(g1,g2))
    ## vector of probabilities for gene g1 to be altered in each sample
    p1 = unlist(lapply(1:length(mutPerSamp), hyperInSamp, g=g1, mutPerSamp=mutPerSamp, nmut=nmut))
    ## vector of probabilities for gene g2 to be altered in each sample
    p2 = unlist(lapply(1:length(mutPerSamp), hyperInSamp, g=g2, mutPerSamp=mutPerSamp, nmut=nmut))
    ## exp12 is a vector with the probability of the g1 and g2 co-occurring in each sample
    exp12 <- p1 * p2
    obs <- length(intersect(which(clonMat[g1,] != 0), which(clonMat[g2,] != 0)))
    rdraws <- randDistrib(exp12, nsim=nsim, mc.cores=15)
    if (obs == 0) {
      if (mean(rdraws) > 0.5) {
        obs=0.5
      } else {
        obs=mean(rdraws)
      }
    }
    es = (obs - mean(rdraws)) / sd(rdraws) ## basically a Z-score
    drawless <- length(which(rdraws < obs))/length(rdraws)

    obsL[[s]][j] <- obs
    unionL[[s]][j] <- length(union(which(clonMat[g1,] != 0), which(clonMat[g2,] != 0))) / ncol(clonMat)
    ratL[[s]][j] = obs / mean(rdraws)
    esL[[s]][j] <- es
    drawp <- ifelse(obs < median(rdraws), length(which(rdraws <= obs))/length(rdraws), length(which(rdraws >= obs))/length(rdraws))
    ## If not observed, set the probability to 1/2 observation in order to avoir error-inducing 0s.
    if (drawp == 0) {
      drawp = 0.5/nsim
    }
    drawpL[[s]][j] <- drawp
    rdrawL[[s]][[j]] <- rdraws
    
    pow2pL[[s]][[g1]][g2]= pwr.2p2n.test(h=es, n1=ncol(clonMat), n2=n1, sig.level=drawp)$power
    pow2pL[[s]][[g2]][g1]= pwr.2p2n.test(h=es, n1=ncol(clonMat), n2=n2, sig.level=alpha)$power
  }
  gc()
}
save(rdrawL, file="rdrawL.RData")
save(topL, file="topL.RData")

## Plot histograms of expect co-occurrences from random draws and actual number of co-occurrences
## for all pairwise interactions in all tumour types
pdf("plots/co-occurrence_expectations.pdf")
for (i in 1:length(sets)) {
  s=sets[i]
  print(s)
  combs <- combL[[s]]
  for (j in 1:length(combs)) {
    tmp <- unlist(strsplit(combs[[j]], split=" - "))
    g1=tmp[1]
    g2=tmp[2]

    rdraws <- rdrawL[[s]][[j]]
    es = esL[[s]][j]
    obs = obsL[[s]][j]
    hist(rdraws, xlim=c(min(c(rdraws,obs)), max(c(rdraws,obs))), main=sprintf("%s: %s - %s (h=%.1f)", s, g1, g2, es),
         xlab="Expected number of co-occurrences")
    abline(v=obs, col="red")
  }
}
dev.off()

library(gplots)
pdf("plots/epistatic_factor_matrices.pdf")
allEpi <- c()
for (s in names(topL)) {
  genes <- topL[[s]]
  combs <- combn(length(genes),2)
  clonMat <- read.table(paste("clonMat_", s, ".txt", sep=""), sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
  fMat <- matrix(NA, ncol=length(genes), nrow=length(genes))
  colnames(fMat) <- rownames(fMat) <- genes
  ## Uncomment to recompute epistatic factor matrices
  ## Currently using the ones from the article. Because random draws are involved, updating these tables can induce variability
  ## fMatZ = fMat
  ## for (i in 1:ncol(combs)) {
  ##   g1 = genes[combs[1,i]]
  ##   g2 = genes[combs[2,i]]
  ##   allEpi[paste(s,g1,g2,sep=":")] <- 1 - ((1 - ratL[[s]][i])) * unionL[[s]][i]
    
  ##   fMat[g1,g2] <- 1 - ((1 - ratL[[s]][i])) * unionL[[s]][i]
  ##   fMat[g2,g1] <- 1 - ((1 - ratL[[s]][i])) * unionL[[s]][i]
  ##   fMatZ[g1,g2] <- esL[[s]][i] * pow2pL[[s]][[g1]][g2]
  ##   fMatZ[g2,g1] <- esL[[s]][i] * pow2pL[[s]][[g2]][g1]
  ##   write.table(fMat, file=paste("epiFactMat_", s, ".txt", sep=""), sep="\t", col.names=NA, quote=F)
  ## }
  ## Comment to recompute epiFactor matrices
  fMat <- read.table(paste("epiFactMat_", s, ".txt", sep=""), sep="\t", header=T, row.names=1, check.names=F)
  heatmap.2(log2(fMat), col=bluered(100), trace="none", dendrogram="none", Colv=NA, Rowv=NA, main=s, mar=c(8,8))
}
dev.off()
