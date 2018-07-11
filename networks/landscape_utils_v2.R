########################################
####           FUNCTIONS            ####
########################################

## Find an optimal set of genes in a matrix of driver genes.
## Simple heuristic that iteratively selects the most frequent gene in the samples not covered yet.
getPanel <- function(panelSize, presMat) {
  panelSize <- min(panelSize, nrow(presMat))
  totPres <- apply(presMat, 1, sum)
  panel <- c()
  subMat <- presMat
  while(length(panel) < panelSize & (!0 %in% dim(subMat))) {
    subPres <- apply(subMat, 1, sum)
    cand <- names(subPres)[which(subPres == max(subPres, na.rm=T))]
    if (length(cand) > 1) {
      cand <- names(totPres[cand])[order(totPres[cand], decreasing=T)[1]]
    }
    panel <- c(panel, cand)
    subMat <- subMat[-which(rownames(subMat)==cand), -which(subMat[cand,] == 1), drop=F]
  }
  ## If all Samples are covered with fewer altereations than the desired panel size, fill with the most frequent ones.
  if (length(panel) < panelSize & nrow(subMat) > 0) {
    sizediff <- min(panelSize - length(panel), nrow(subMat))
    panel <- c(panel, names(sort(totPres[rownames(subMat)], decreasing=T))[1:sizediff])
  }
  return(panel)
}

coOccurNratio <- function(g1, g2, presMat, gprob, alpha=0.05, nrand=1000, pow=0.8, maxn=10000) {
  if (g1==g2) {
    return(1)
  }
  n <- ncol(presMat)
  n1 <- sum(presMat[g1,])
  n2 <- sum(presMat[g2,])
  n12 <- length(which(apply(presMat[c(g1,g2),],2,sum) == 2))
  p1 = n2/n
  p12 = n12/n1

  ##print(c(g1,g2))
  h = ES.h(p1, p12) ##(2 * asin(sqrt(p1))) - (2 * asin(sqrt(p12)))
  idealn <- tryCatch(pwr.2p2n.test(h=h, power=pow, n1=n, sig.level=alpha)$n2, error=function(err){n2*maxn})
  return(max(min(n2/idealn, 1), 1/maxn))
}

## Ratio between observed and expected
coOccurR <- function(g1, g2, presMat, gprob, pbefore, totMut=NULL) {
  if (is.null(totMut)) {
    totMut <- apply(presMat, 2, sum)
  }
  expBinA <- unlist(lapply((totMut-1)[which(presMat[g1,] == 1)], function(n, p){(1 -((1 - p) ^ n))}, p=gprob[g2]))
  expAinB <- unlist(lapply((totMut-1)[which(presMat[g2,] == 1)], function(n, p){(1 -((1 - p) ^ n))}, p=gprob[g1]))
  expAB <- (sum(expBinA) * pbefore[g1,g2]) + (sum(expAinB) * pbefore[g2,g1])
  obsAB <- length(which(apply(presMat[c(g1,g2),], 2, sum) == 2))
  if (obsAB == 0) 
    obsAB <- min(0.5, expAB/2)
  ##naUb <- length(which(apply(presMat[c(g1,g2),], 2, sum) > 0))
  return(obsAB/expAB)
  ##return(exp(zab * (naUb / ncol(presMat))))
}

## Binomial-based power for the number of co-occurences (based on detecting a doubled probability)
coOccurPow <- function(g1, g2, presMat, gprob, alpha=0.05, nrand=1000, r=2) {
  n <- ncol(presMat)
  n1 <- sum(presMat[g1,])
  n2 <- sum(presMat[g2,])
  n12 <- length(which(apply(presMat[c(g1,g2),],2,sum) == 2))
  p1 = n2/n
  p12 = n12 / n1

  h = ES.h(p1, ifelse(p12 >= p1, min(r*p1,1), p1/r)) ##(2 * asin(sqrt(p1))) - (2 * asin(sqrt(p12)))
  return(tryCatch(pwr.2p2n.test(h=h, n1=n, n2=n1, sig.level=alpha)$power, error=function(err){0}))

  ## rand12 <- unlist(lapply(1:nrand, function(x){length(which(sample(presMat[g2,], size=n1, replace=T)==1))}))
  ## randp <- unlist(lapply(rand12, function(x){phyper(x, n2, n-n2, n1)}))

  ##   totMut <- apply(presMat, 2, sum)
  ## }
  ## expAB <- unlist(lapply(totMut, function(n, p1, p2){(1 -((1 - p1) ^ n)) * (1 -((1 - p2) ^ n))}, p1=gprob[g1], p2=gprob[g2]))
  ## return(tryCatch(pwr.2p2n.test(h=(sum(expAB) - n12), n1=n1, n2=n2, sig.level=alpha)$power, error=function(err){0}))
}

## Epistatis factor
coOccurZ <- function(g1, g2, presMat, gprob, totMut=NULL) {
  if (is.null(totMut)) {
    totMut <- apply(presMat, 2, sum)
  }
  expAB <- unlist(lapply(totMut, function(n, p1, p2){(1 -((1 - p1) ^ n)) * (1 -((1 - p2) ^ n))}, p1=gprob[g1], p2=gprob[g2]))
  obsAB <- length(which(apply(presMat[c(g1,g2),], 2, sum) == 2)) / ncol(presMat)
  ##naUb <- length(which(apply(presMat[c(g1,g2),], 2, sum) > 0))
  return((obsAB - mean(expAB)) / sd(expAB))
  ##return(exp(zab * (naUb / ncol(presMat))))
}

getStartFit <- function(g, presMat, gprob, totMut=NULL) {
  return(sum(presMat[g,]) / ncol(presMat))
  ## if (is.null(totMut)) {
  ##   totMut <- apply(presMat, 2, sum)
  ## }
  ## expg <- unlist(lapply(totMut, function(n, p){1 -((1 - p) ^ n)}, p=gprob[g]))
  ## obsg <- sum(presMat[g,]) / ncol(presMat)
  ## zg <- (obsg - mean(expg)) / sd(expg)
  ## return(mean(expg) + (zg * mean(expg)))
}

getFit <- function(subpanel, fit, ints) {
  return(sum(unlist(lapply(subpanel, function(g){(prod(ints[g,setdiff(subpanel,g)]) * fit[g])}))))
}

getFitStr <- function(panelStr, fit, ints) {
  subpanel=unlist(strsplit(panelStr, ":"))
  return(getFit(subpanel, fit, ints))
}

getCombFit <- function(allComb, panel, fit, ints) {
  combFit <- unlist(apply(allComb, 2, function(gidx){getFit(panel[gidx], fit, ints)}))
  names(combFit) <- apply(allComb, 2, function(x){paste(panel[x],collapse=":")})
  return(combFit)
}

getNsamp <- function(subpanel, presMat, gprob, alpha=0.05, power=0.8, r=2, totMut=NULL) {
  ##print(subpanel)
  if (is.null(totMut)) {
    totMut <- apply(presMat, 2, sum)
  }
  expComb <- unlist(lapply(totMut, function(n, p){prod(unlist(lapply(p, function(x){(1-((1 - x) ^ n))})))}, p=gprob[subpanel]))
  prob <- sum(expComb) / ncol(presMat)
  es <- ES.h(prob,prob*r)
  return(tryCatch(pwr.p.test(h=es, sig.level=alpha, power=power)$n), error=function(err){NA})
}

combNsample <- function(allComb, panel, presMat, gprob, alpha=0.05, power=0.8, r=2, totMut=NULL) {
  print(dim(allComb))
  combNsamp <- unlist(apply(allComb, 2, function(gidx){getNsamp(panel[gidx], presMat=presMat, gprob=gprob, alpha=alpha, power=power, r=r, totMut=totMut)}))
  names(combNsamp) <- apply(allComb, 2, function(x){paste(panel[x],collapse=":")})
  return(combNsamp)
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
