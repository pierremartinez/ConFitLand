## To be run in the tcga/ directory

load("./topL.RData")

jacSimMat <- function(m) {
  interidx = which(apply(m, 1, sum) == 2)
  interl = length(interidx)
  ##noninter = setdiff(1:nrow(m), interidx)
  a = sum(m[,1])
  b = sum(m[,2])
  if ((a+b) == 0)
    return(1)
  interl / (a + b - interl)
}

jacSimIdx <- function(idx,idxComb,m) {
  jacSimMat(m[,c(idxComb[idx,1],idxComb[idx,2])])
}

nsim=10
library(parallel)
realJacL <- simJacL <- list()
for (s in names(topL)) {
  print(s)
  clonMat <- read.table(paste("./clonMat_", s, ".txt", sep=""), sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
  ## Only considered whether the gene is mutated or not in a 0/1 matrix
  clonMat01 <- clonMat
  clonMat01[clonMat01 > 1] <- 1
  idxComb <- matrix(nr=(ncol(clonMat01) * (ncol(clonMat01) - 1) / 2), nc=2)
  idx <- 1
  for (i in 1:(ncol(clonMat01)-1)) {
    for (j in (i+1):ncol(clonMat01)) {
      idxComb[idx,] <- c(i,j)
      idx=idx+1
    }
  }
  ## Calculate Jaccard Indices
  realJac <- unlist(mclapply(1:nrow(idxComb), jacSimIdx, idxComb=idxComb, m=clonMat01, mc.cores=15))
  realJacL[[s]] <- realJac

  ## Build randomised matrices with similar distributions of clonal alterations per patient
  simJacL[[s]] <- list()
  for (sim in 1:nsim) {
    print(sim)
    matsim <- matrix(0,nr=nrow(clonMat01), nc=ncol(clonMat01))
    for (i in 1:ncol(matsim)) {
      matsim[sample(1:nrow(clonMat01), sum(clonMat01[,i])),i] <- 1
    }
    ## Calculate Jaccard Indices
    simJacL[[s]][[sim]] <- unlist(mclapply(1:nrow(idxComb), jacSimIdx, idxComb=idxComb, m=matsim, mc.cores=15))
  }
}
save(realJacL, file="realJacL.RData")
save(simJacL, file="simJacL.RData")

ratio95 <- function(s, realL, simL) {
  simpool <- simL[[s]]
  sim95 <- quantile(simpool, .95, na.rm=T)
  realL[[s]] / sim95
}

simpoolL <- lapply(simJacL, unlist)
wtp <- unlist(lapply(sets, function(x){wilcox.test(realJacL[[x]], simpoolL[[x]])$p.value}))
realRatio95 <- lapply(names(realJacL), ratio95, realL=realJacL, simL=simpoolL)
perc95 <- unlist(lapply(realRatio95, function(x){length(which(x>1)) / length(x)}))
names(perc95) <- names(realRatio95) <- names(realJacL)

rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")
pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- names(socPowL)
colSet <- rainbow10[1:length(socPowL)]
names(colSet) <- names(socPowL)

library(beanplot)
pdf("plots/repeatability_plots.pdf")
##beanplots
beanplot(realRatio95, cutmin=0, what=c(0,1,1,0), col=rainbow11[c(8,3,4,1)], show.names=T, las=2, ylab="Normalised Jaccard Index")
abline(h=1, lty=2)
text(unlist(lapply(perc95*100, sprintf, fmt="%.0f%%")), x=seq(0.66, length(sets)-.34, 1), y=max(unlist(realRatio95), na.rm=T)/3,
     col=colSet[sets])

beanplot(realJacL, cutmin=0, what=c(0,1,1,0), col=c(rainbow11[c(8,3,4)], "black"), xlim=c(0.5, length(realJacL)*3-.5), bw=.001,
         ylim=c(0, max(c(unlist(realJacL), unlist(simJacL)), na.rm=T)), at=seq(1, length(realJacL)*3-2, 3), show.names=F,
         ylab="Jaccard Index", maxwidth=1.25)
beanplot(simpoolL, cutmin=0, what=c(0,1,1,0), col=c(rainbow11[c(8,3,4)], "black"), bw=.001, at=seq(2, length(simpoolL)*3-1, 3), add=T,
         show.names=F, maxwidth=1.25)
axis(1, at=seq(1.5, length(realJacL)*3-1.5, 3), label=sets, las=2)

##boxplots
boxplot(realRatio95, outline=F, col=rainbow10, las=2, ylab="Normalised Jaccard Index")
abline(h=1, lty=1)
text(unlist(lapply(perc95*100, sprintf, fmt="%.0f%%")), x=seq(0.66, length(sets)-.34, 1), y=7.25,
     col=colSet[sets])

boxplot(realJacL, cutmin=0, what=c(0,1,1,0), col=rainbow10, xlim=c(0.5, length(realJacL)*3-.5), bw=.001,
         at=seq(1, length(realJacL)*3-2, 3), xaxt="n", ylab="Jaccard Index", maxwidth=1.25, outline=F)
boxplot(simpoolL, cutmin=0, what=c(0,1,1,0), col=rainbow10, bw=.001, at=seq(2, length(simpoolL)*3-1, 3), add=T,
         show.names=F, maxwidth=1.25, outline=F)
axis(1, at=seq(1.5, length(realJacL)*3-1.5, 3), label=sets, las=2)
dev.off()

unlist(lapply(realRatio95, mean))
