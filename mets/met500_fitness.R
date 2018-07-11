library(gdata)
files <- dir("../tcga", pattern="clonMat")
sets <- toupper(gsub("^clonMat_(\\w+)\\.txt$", "\\1", files))

#### Get IntOGen drivers
intomut <- read.table("../misc/Mutational_drivers_per_tumor_type.tsv", sep="\t", header=T, stringsAsFactors=F, check.names=F)
intocna <- read.table("../misc/CNA_drivers_per_tumor_type.tsv", sep="\t", header=T, stringsAsFactors=F, check.names=F)
equivSet <- c("BLCA", "BRCA", "COREAD", "GBM", "HNSC", "RCCC", "LUAD", "LUSC", "CM")
names(equivSet) <- sets
drivL <- list()
for (s in sets) {
  tmp <- paste(intocna$geneHGNCsymbol, intocna$GISTIC_CNA, sep="_")[which(intocna$Tumor_type_GISTIC == equivSet[s])]
  tmp <- gsub("_A$", "_gain", gsub("_D$", "_loss", tmp))
  drivL[[s]] <- c(intomut$geneHGNCsymbol[which(intomut$Tumor_type == equivSet[s])], tmp)
}

patDat <- read.xls("sample_info.xlsx", stringsAsFactors=F)

## Attributing a tumour type to each reported case
cantype <- c("COAD",
"SKCM",
NA,
NA,
NA,
NA,
NA,
"BLCA",
NA,
NA,
NA,
NA,
NA,
NA,
"BRCA",
"BRCA",
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
"BLCA",
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
"KIRC",
NA,
"COAD",
NA,
NA,
NA,
NA,
NA,
NA,
"BLCA",
NA,
NA,
NA,
"SKCM",
"LUAD",
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
"LUSC",
"HNSC",
NA,
NA,
NA,
NA,
"KIRC",
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
"BRCA",
NA,
"HNSC",
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
"HNSC",
NA,
"BRCA",
NA,
"GBM",
NA,
NA,
"Myoepithelial Carcinoma",
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
"BRCA",
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
NA,
"HNSC")

names(cantype) <- c("Colon Adenocarcinoma",
"Cutaneous Melanoma",
"Solitary Fibrous Tumor",
"Extrahepatic Cholangiocarcinoma",
"Prostate Adenocarcinoma",
"Dedifferentiated Liposarcoma",
"Prostate Neuroendocrine Carcinoma",
"Bladder Urothelial Carcinoma",
"Thymic Neuroendocrine Tumor",
"Extraskeletal Myxoid Chondrosarcoma",
"Osteoblastic Osteosarcoma",
"Embryonal Rhabdomyosarcoma",
"Undifferentiated Pleomorphic Sarcoma",
"Adrenocortical Carcinoma",
"Breast Invasive Lobular Carcinoma",
"Breast Invasive Ductal Carcinoma",
"Anal Squamous Cell Carcinoma",
"Atypical Lung Carcinoid",
"Gastrointestinal Stromal Tumor",
"Desmoplastic Small-Round-Cell Tumor",
"Intrahepatic Cholangiocarcinoma",
"Salivary Clear Cell Carcinoma",
"Adenoid Cystic Carcinoma",
"Pancreatic Adenocarcinoma",
"Penile Squamous Cell Carcinoma",
"Ewing Sarcoma",
"Bladder Squamous Cell Carcinoma",
"Serous Ovarian Cancer",
"Anaplastic Thyroid Cancer",
"Cutaneous Squamous Cell Carcinoma",
"Follicular Dendritic Cell Sarcoma",
"Pancreatic Neuroendocrine Tumor",
"High-Grade Serous Ovarian Cancer",
"Thymoma",
"Ameloblastoma",
"Renal Clear Cell Carcinoma",
"Esophageal Adenocarcinoma",
"Rectal Adenocarcinoma",
"Poorly Differentiated Carcinoma, NOS",
"Dermatofibrosarcoma Protuberans",
"Poorly Differentiated Thyroid Cancer",
"Sarcomatoid Carcinoma",
"Undifferentiated Malignant Neoplasm",
"Solid Pseudopapillary Neoplasm of the Pancreas",
"Small Cell Bladder Cancer",
"Mucinous Adenocarcinoma",
"Leiomyosarcoma",
"Thymic Squamous Carcinoma",
"Head and Neck Mucosal Melanoma",
"Lung Adenocarcinoma",
"Merkel Cell Carcinoma",
"Pancreatobiliary Ampullary Carcinoma",
"Hepatocellular Carcinoma",
"Poorly Differentiated Non-Small Cell Lung Cancer",
"Mucinous Cystic Neoplasm",
"Synovial Sarcoma",
"Fibrolamellar Carcinoma",
"Uterine Leiomyosarcoma",
"Squamous Cell Carcinoma, NOS",
"Clear Cell Carcinoma",
"Alveolar Soft Part Sarcoma",
"Teratoma with Malignant Transformation",
"Adenocarcinoma, NOS",
"Lung Squamous Cell Carcinoma",
"Oral Cavity Squamous Cell Carcinoma",
"Low-Grade Serous Ovarian Cancer",
"Urachal Adenocarcinoma",
"Undifferentiated Pleomorphic Sarcoma/Malignant F",
"Clear Cell Carcinoma of the Lung",
"Renal Cell Carcinoma",
"Urethral Urothelial Carcinoma",
"Angiosarcoma",
"Skin Adnexal Carcinoma",
"Myxoid/Round-Cell Liposarcoma",
"Stomach Adenocarcinoma",
"Mullerian Poorly Differentiated Carcinoma",
"Chromophobe Renal Cell Carcinoma",
"Small Cell Carcinoma of Unknown Primary",
"Perivascular Epithelioid Cell Tumor",
"Neuroendocrine Carcinoma of the Stomach",
"Inflammatory Breast Cancer",
"Yolk Sac Tumor",
"Oropharynx Squamous Cell Carcinoma",
"Soft Tissue Myoepithelial Carcinoma",
"Myxoid Chondrosarcoma",
"Embryonal Carcinoma",
"Signet Ring Cell Carcinoma of the Stomach",
"Salivary Duct Carcinoma",
"Granular Cell Tumor",
"Mullerian Adenocarcinoma ",
"High-Grade Sarcoma",
"Esophageal Squamous Cell Carcinoma",
"Adamantinoma",
"Metaplastic Breast Cancer",
"Ovarian Adenocarcinoma",
"Glioblastoma Multiforme",
"Myxofibrosarcoma",
"Rhabdomyosarcoma",
"Myoepithelial Carcinoma",
"Sex Cord Stromal Tumor",
"Head and Neck Neuroendocrine Carcinoma",
"Adenosarcoma",
"Neuroendocrine Tumor, NOS",
"Collecting Duct Renal Cell Carcinoma",
"Osteosarcoma",
"Salivary Carcinoma",
"Salivary Carcinoma, Other",
"Cervical Neuroendocrine Tumor",
"Breast Invasive Carcinoma, NOS",
"Uterine Endometrioid Carcinoma",
"Undifferentiated High-Grade Sarcoma",
"Small Cell Lung Cancer",
"Sebaceous Carcinoma",
"Seminal Vesicle Carcinoma",
"Atypical Meningioma",
"Extramammary Paget Disease",
"Prostate Carcinosarcoma",
"Clear Cell Ovarian Cancer",
"Pleural Mesothelioma",
"Dedifferentiated Chondrosarcoma",
"Ependymoma",
"Sinonasal Squamous Cell Carcinoma")

## Attributing a type to each metastatic location. LUAD as default for lung, change to LUSC if desired.
lungtype="LUAD"
biopsy <- c(NA, "SKCM", NA, NA, NA, NA,
            lungtype, NA, NA, NA, NA, NA,
            "BRCA", NA, NA, NA, NA, "HNSC", NA, NA, NA,
            "COAD", "HNSC", NA, NA, NA, NA,
            "HNSC", "GBM", "BLCA", NA, NA, NA, 
            NA, NA, NA, NA, NA, NA, 
            NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA,
            "SKCM", NA, NA, NA, NA, "HNSC",
            NA, NA, NA, NA, NA, NA,
            "HNSC", NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA,
            NA, NA)
names(biopsy) <- c("Liver", "Skin", "Dura", "Shoulder", "Subq Nodule", "Lymph Node",
                   "Lung", "Chest Wall", "Subq Abdom. Mass", "Forearm Mass", "Retroperitoneum", "Pancreas",
                   "Breast", "Bone", "Gluteus Muscle","Thyroid", "Adrenal", "Neck", "Vertebra", "Esophagus", "Buttock ",
                   "Colon", "Mandible", "Periaortic Mass", "Abdom. Wall", "Thoracic Epidural", "Paraspinal Mass",
                   "Cheek", "Brain", "Bladder", "Epidural", "Rectus Muscle", "Pelvis", 
                   "Pelvic Mass", "Peritoneum", "Pleura", "Abdomen", "Gluteal Mass", "Subq Nodule (Arm)", 
                   "Peritoneal Fluid", "Subcut. Abdomen", "Subcut. Periumbilical","Acetabulum", "Lung Nodule", "Perirectal Mass",
                   "Thigh", "Chest Wall Nodule", "Parotid", "Omentum", "Hilum", "Abdom. Nodule",
                   "Skin Punch", "Sternal Mass", "Infraclavicular Mass", "Abdom. Mass", "Subq", "Neck Mass",
                   "Subq Chest Wall", "Shoulder Mass", "Pleural Mass", "Back", "Thigh Mass", "Scapular Mass",
                   "Submandible", "Mesenteric Mass", "Peritoneal Nodule", "Psoas Muscle", "Thoracic Fluid", "Hip Lesion",
                   "Axillary Mass", "Omental Mass", "Cervix", "Spinal Mass", "Lacrimal Gland", "Prostate",
                   "Arm", "Sphenoid Sinus")

## Reformat patient data
patDat$Met <- biopsy[patDat$BiopsySite]
patDat$Prim <- cantype[patDat$CancerType]
save(patDat, file=paste("patDat_",lungtype,".RData", sep=""))

## Load sequencing info (sample purity estimation)
seqDat <- read.xls("sequencing_info.xlsx", stringsAsFactors=F)
purity <- as.numeric(gsub("\\%", "", gsub("<30\\%", "25%", seqDat$Sequencing))) / 100
names(purity) <- seqDat$Case

## Load and reformat mutation data
mutDat <- read.table("somatic_v4.csv", header=T, sep=",", stringsAsFactors=F, check.names=F)
mutDat$Pat <- gsub("^(\\w+)\\.SI.+$", "\\1", mutDat$Pipeline_ID)
mutDat$VAF <- mutDat$Var_T / mutDat$Tot_T
mutDat$Chr[mutDat$Chr == "X"] <- 23
mutDat$Chr[mutDat$Chr == "Y"] <- 24
mutDat$Chr = as.numeric(mutDat$Chr)
mutDat$Effect <- gsub("^([^;]+);.+$", "\\1", mutDat$Effect, perl=T)
mutDat$Effect[mutDat$Effect == "Disruptive In-Frame Deletion"] <- "Missense_Mutation"
mutDat$Effect[mutDat$Effect == "Disruptive In-Frame Insertion"] <- "Missense_Mutation"
mutDat$Effect[mutDat$Effect == "Exon Loss"] <- "Frame_Shift_Del"
mutDat$Effect[mutDat$Effect == "Frameshift"] <- "Frame_Shift_Del"
mutDat$Effect[mutDat$Effect == "In-Frame Deletion"] <- "Missense_Mutation"
mutDat$Effect[mutDat$Effect == "In-Frame Insertion"] <- "Missense_Mutation"
mutDat$Effect[mutDat$Effect == "Missense"] <- "Missense_Mutation"
mutDat$Effect[mutDat$Effect == "Splice Acceptor"] <- "Splice_site"
mutDat$Effect[mutDat$Effect == "Splice Donor"] <- "Splice_site"
mutDat$Effect[mutDat$Effect == "Startloss"] <- "Nonsense_Mutation"
mutDat$Effect[mutDat$Effect == "Stopgain"] <- "Nonsense_Mutation"
mutDat$Effect[mutDat$Effect == "Stoploss"] <- "Nonstop_Mutation"

## Load and reformat CNA data
cnvDat <- read.table("cnv_v4.csv", header=T, sep=",", stringsAsFactors=F, check.names=F)
cnvDat$Pat <- gsub("^(\\w+)\\.SI.+$", "\\1", cnvDat$Pipeline_ID)
ab <- gsub("^.*Genotype *", "", cnvDat$Classification)
nanb <- matrix(nc=2, unlist(lapply(ab, function(x) {
  if (x == "")
    return(c(1,1))
  l <- unlist(strsplit(x, split="/"))
  return(c(length(which(l == "A")), length(which(l == "B"))))
})), byrow=T)
pbidx <- which(apply(nanb, 1, sum) == 0 & cnvDat$Copy_Number > 0)
nanb[pbidx,2] <- cnvDat$Copy_Number[pbidx]
cnvDat$nA <- nanb[,1]
cnvDat$nB <- nanb[,2]
cnvDat$Chr[cnvDat$Chr == "X"] <- 23
cnvDat$Chr[cnvDat$Chr == "Y"] <- 24
cnvDat$Chr = as.numeric(cnvDat$Chr)

## Calculate the ploidy of each sample
ploidy <- unlist(lapply(unique(cnvDat$Pat), function(p){
  subp <- cnvDat[which(cnvDat$Pat == p),]
  sum(subp$Copy_Number* subp$Targeted_Exons) / sum(subp$Targeted_Exons)
}))
names(ploidy) <- unique(cnvDat$Pat)

## Which patients have 1 or 2 
validPats <- intersect(patDat$ID[which(patDat$Prim %in% sets | patDat$Met %in% sets)],
                       intersect(unique(cnvDat$Pat), unique(mutDat$Pat)))
doublePats <- intersect(patDat$ID[which(patDat$Prim %in% sets & patDat$Met %in% sets)],
                        intersect(unique(cnvDat$Pat), unique(mutDat$Pat)))
subMut <- mutDat[which(mutDat$Pat %in% validPats & mutDat$Chr %in% as.character(1:22)),]
subCnv <- cnvDat[which(cnvDat$Pat %in% validPats & cnvDat$Chr %in% as.character(1:22)),]

segments <- data.frame(SampleID = subCnv$Pat,
                       Chr = as.numeric(subCnv$Chr),
                       Start = subCnv$Start,
                       End = as.numeric(subCnv$End),
                       nProbes = as.numeric(subCnv$Targeted_Exons),
                       cn = as.numeric(subCnv$Copy_Number),
                       nA = as.numeric(subCnv$nA),
                       nB = as.numeric(subCnv$nB),
                       Ploidy = as.numeric(ploidy[subCnv$Pat]),
                       `Aberrant Cell Fraction` = as.numeric(purity[subCnv$Pat]),
                       ChipNames = subCnv$Pipeline_ID,
                       NormalFile = "None",
                       TumorTissue = "Metastasis",
                       NormalTissue = "Blood",
                       check.names=F, stringsAsFactors=F)
met500cn <- list(segments=segments[which(segments$cn > 0),])
save(met500cn, file=paste("met500_cn.RData", sep=""))

mutTab <- data.frame(Patient = subMut$Pat,
                     Chr = as.numeric(subMut$Chr) ,
                     Start_position = subMut$Pos ,
                     End_position = subMut$Pos + nchar(subMut[,"substr(Ref,1,100)"]) - 1,
                     Reference = subMut[,"substr(Ref,1,100)"],
                     Alternate = subMut[,"substr(Alt,1,100)"],
                     Variant_freq = subMut$Var_T,
                     Ref_freq = subMut$Tot_T - subMut$Var_T,
                     Hugo_Symbol = subMut$Gene,
                     Variant_Classification = subMut$Effect,
                     Protein_Change = subMut[,"substr(Protein_Change,1,100)"])

## library(parallel)
## mutInCNV <- unlist(mclapply(1:nrow(mutTab), function(i) {
##   x <- mutTab[i,]
##   subp <- segments[which(segments$SampleID == x$Pat & segments$Chr == x$Chr),]
##   return(length(which(subp$Start < x$End & subp$End > x$Start)))
## }, mc.cores=8))
## mutTab=mutTab[which(mutInCNV > 0),]
write.table(mutTab, file="met500_mutation_table.txt", sep="\t", quote=F, row.names=F)

## Calculate CCF
outdir <- "met500_ccf"
if (!dir.exists(outdir))
  dir.create(outdir)
## Package available from http://stm.sciencemag.org/content/7/283/283ra54.short
library("EstimateClonality")
for (p in validPats) {
  tmp <- clonality.estimation(mutation.table.loc="met500_mutation_table.txt"
                             ,seg.mat.loc=paste("met500_cn.RData", sep="")
                             ,data.type=outdir
                             ,TCGA.barcode=p
                             ,ANALYSIS.DIR="./")
}

id2gene <- as.character(mutTab$Hugo_Symbol)
names(id2gene) <- gsub(" ", "", apply(mutTab[,c("Patient","Chr","Start_position","Reference")], 1, paste, collapse=":"))

## Investigate CN events
cne <- unlist(drivL)
cne <- unique(cne[grep("_", cne)])
cneg <- unique(unlist(lapply(cne, function(x){unlist(strsplit(x, "_"))[1]})))
library(biomaRt)
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)
ensres <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                filter="hgnc_symbol", values=cneg, mart=ensembl)

library(gdata)
glDat <- read.xls("nature23306-s3.xlsx", skip=1)
equivCn <- c("amplification", "deletion")
names(equivCn) <- c("gain", "loss")
cnMat <- matrix(0, nr=length(cne), nc=length(validPats))
colnames(cnMat) <- validPats
rownames(cnMat) <- cne
for (cn in cne) {
  l <- unlist(strsplit(cn, split="_"))
  st <- as.numeric(ensres[which(ensres$hgnc_symbol == l[1]),"start_position"])
  e <- as.numeric(ensres[which(ensres$hgnc_symbol == l[1]),"end_position"])
  chr <- ensres[which(ensres$hgnc_symbol == l[1]),"chromosome_name"]
  if (chr == "X")
    chr=23
  if (chr == "Y")
    chr=24
  chr=as.numeric(chr)
  subsegs <- segments[which(segments$Chr == chr),]
  subsegs <- subsegs[which(subsegs$Start < e & subsegs$End > st),]
  if (l[2] == "gain") {
    subsegs <- subsegs[which(subsegs$cn > (subsegs$Ploidy + 0.6)),]
  } else {
    subsegs <- subsegs[which(subsegs$cn < (subsegs$Ploidy - 0.6) | subsegs$nA == 0 | subsegs$nB == 0),]
  }
  cnMat[cn,unique(subsegs$SampleID)] <- 1
  ## subgl <- glDat[which(glDat$gene == l[1] & glDat$variant_class == equivCn[l[2]]),]
  ## cnMat[cn,intersect(subgl$sample, validPats)] <- 1
}

## Build clonMats
d <- dir(outdir)
allMut <- c()
for (s in d) {
  sdat <- read.table(paste(outdir, "/", s, "/", s, ".earlylate.tsv", sep=""), sep="\t", check.names=F, header=T, stringsAsFactors=F)
  allMut <- rbind(allMut, sdat)
}
allMut <- cbind(allMut, as.character(id2gene[allMut[,"mutation_id"]]))
colnames(allMut)[ncol(allMut)] <- "gene"
write.table(allMut, file="met500_allCCF.txt", sep="\t", row.names=F, quote=F)

clonMut <- allMut[which(allMut[,"absolute.ccf.0.95"] == 1),]
allg <- sort(unique(clonMut[,"gene"]))
clonMat <- matrix(0, nr=length(allg), nc=length(validPats))
rownames(clonMat) <- allg
colnames(clonMat) <- validPats
for (i in 1:nrow(clonMut)) {
  clonMat[as.character(clonMut[i,"gene"]), as.character(clonMut[i,"patient"])] <-
    clonMat[as.character(clonMut[i,"gene"]), as.character(clonMut[i,"patient"])] + 1
}
clonMat <- rbind(clonMat, cnMat)
write.table(clonMat, file="met500_clonMat.txt", sep="\t", col.names=NA, quote=F)

################################################
####           Fitness Calculation          ####
################################################

lungtype="LUAD"

load(paste("../patDat_",lungtype,".RData", sep=""))
clonMat <- read.table(paste("../met500_clonMat.txt",sep=""), header=T, sep="\t", check.names=F, stringsAsFactors=F, row.names=1)
validPats <- colnames(clonMat)
prim <- patDat$Prim
met <- patDat$Met
names(prim) <- names(met) <- patDat$ID

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

## Data on select CN events
load("../tcga/cnSelectL.RData")

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

maxland = 50
medFit <- unlist(lapply(bestFitL, median, na.rm=T))
names(medFit) <- names(bestFitL)
fitMat <- matrix(nr=length(validPats), nc=4+length(sets))
rownames(fitMat) <- validPats
colnames(fitMat) <- c("Prim", "Met", "PrimFit", "MetFit", sets)
fitMat[,"Prim"] <- prim[validPats]
fitMat[,"Met"] <- met[validPats]
library(parallel)
## Compute fitness in all landscapes (i.e. for each TCGA set)
for (s in sets) {
  print(s)
  intnf = bestCombs[s,1]
  epinf = bestCombs[s,2]
  genes <- topL[[s]][1:min(length(topL[[s]]), maxland)]
  fMat <- read.table(paste("../tcga/epiFactMat_", s, ".txt", sep=""),
                     sep="\t", header=T, row.names=1, check.names=F)
  fMat <- fMat[genes,]
  sMat <- read.table(paste("../tcga/clonMat_", s, ".txt", sep=""),
                     sep="\t", row.names=1, stringsAsFactors=F, header=T, check.names=F)
  nmut <- apply(sMat, 1, sum)
  selectMut <- c(selectMutL[[s]], cnSelectL[[s]][rownames(sMat)[grep("_", rownames(sMat))]])
  genes <- intersect(genes, rownames(clonMat))
  ssMat <- clonMat[genes,]
  mutPerSamp <- apply(ssMat, 2, sum)
  socRat <- (1 - (socRatL[[s]][genes] - 1) * socPowL[[s]][genes])
  names(socRat) <- genes
  geneComb <- apply(ssMat, 2, function(x){paste(rownames(ssMat)[which(x==1)], collapse=":")})
  sFit <- unlist(mclapply(geneComb, getFitStr, sel=selectMut, socRat=socRat, ints=fMat, sf=intnf, ef=epinf, mc.cores=8))
  sFitN <- sFit / median(medFit[s], na.rm=T)

  fitMat[validPats,s] <- sFitN
}

## Define fitness in primary and metastatic landscapes
for (i in 1:nrow(fitMat)) {
  if (fitMat[i,"Prim"] %in% sets) {
    fitMat[i,"PrimFit"] <- fitMat[i,fitMat[i,"Prim"]]
  }
  if (fitMat[i,"Met"] %in% sets) {
    fitMat[i,"MetFit"] <- fitMat[i,fitMat[i,"Met"]]
  }
}

library(beeswarm)
## Samples with paired primary and metastatic landscapes
bothIdx <- which(!is.na(fitMat[,"PrimFit"]) & !is.na(fitMat[,"MetFit"]))
rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")
pchSet <- c(16,15,17,18,1,0,2,5,3)
names(pchSet) <- names(socPowL)
colSet <- rainbow10[1:length(socPowL)]
names(colSet) <- names(socPowL)

pdf(paste("plots/met500_fitness_",lungtype,".pdf",sep=""))
wt <- wilcox.test(as.numeric(fitMat[bothIdx,"MetFit"]), as.numeric(fitMat[bothIdx,"PrimFit"]), paired=T)
boxplot(as.numeric(fitMat[bothIdx,"PrimFit"]), as.numeric(fitMat[bothIdx,"MetFit"]), ylab="Normalised Fitness",
        sub=sprintf("p=%.3f (paired Wilcoxon test)", wt$p.value), main="Cases with known Prim and Met fitness",
        outline=F, ylim=c(0, max(c(as.numeric(fitMat[bothIdx,"MetFit"]), as.numeric(fitMat[bothIdx,"PrimFit"])))),
        names=c("Primary", "Metastatasis"))
beeswarm(c(as.numeric(fitMat[bothIdx,"PrimFit"]), as.numeric(fitMat[bothIdx,"MetFit"])) ~ rep(1:2, rep(length(bothIdx), 2)),
         add=T, cex=1.5, pwpch=pchSet[c(fitMat[bothIdx,"Prim"], fitMat[bothIdx,"Met"])],
         pwcol=colSet[c(fitMat[bothIdx,"Prim"], fitMat[bothIdx,"Met"])], corral="wrap")
for (i in bothIdx) {
  arrows(1, as.numeric(fitMat[i,"PrimFit"]), 2, as.numeric(fitMat[i,"MetFit"]), code=0)
}
legend("topright", legend=names(bestFitNL), col=colSet[names(bestFitNL)], pch=pchSet[names(bestFitNL)])

wt <- wilcox.test(as.numeric(fitMat[,"PrimFit"]), as.numeric(fitMat[,"MetFit"]))
boxplot(as.numeric(fitMat[,"PrimFit"]), as.numeric(fitMat[,"MetFit"]), sub=sprintf("p=%.3f (Wilcoxon test)", wt$p.value),
        names=c("Primary", "Metastasis"), ylab="Normalised fitness", outline=F,
        ylim=c(0, max(as.numeric(fitMat[,c("PrimFit","MetFit")]), na.rm=T)))
beeswarm(as.numeric(fitMat[,c("PrimFit","MetFit")]) ~ rep(1:2, rep(length(validPats),2)), add=T, cex=1.5, corral="wrap",
         pwpch=pchSet[c(fitMat[,"Prim"], fitMat[,"Met"])], pwcol=colSet[c(fitMat[,"Prim"], fitMat[,"Met"])])
legend("topright", legend=names(bestFitNL), col=colSet[names(bestFitNL)], pch=pchSet[names(bestFitNL)])

dev.off()

## Mean fitness in the other landscapes
mOtherMet <- apply(fitMat, 1, function(x){
  tmps <- setdiff(sets, x[c("Met", "Prim")])
  if (length(tmps) == length(sets))
    return(NA)
  mean(as.numeric(x[tmps]))
})

## Metastases
library(beeswarm)
pdf(paste("plots/mets_VS_otherMets_",lungtype, "_", model, ".pdf", sep=""))
metIdx <- which(!is.na(fitMat[,"MetFit"]))
wt <- wilcox.test(as.numeric(fitMat[metIdx,"MetFit"]), mOtherMet[metIdx], paired=T)
boxplot(as.numeric(fitMat[metIdx,"MetFit"]), mOtherMet[metIdx], ylab="Normalised fitness",
        sub=sprintf("p=%.3f (paired Wilcoxon test)", wt$p.value), main="Metastatic cases",
        outline=F, ylim=c(0, max(c(as.numeric(fitMat[metIdx,"MetFit"]), mOtherMet[metIdx]))),
        names=c("Specific landscape", "Other landscapes (mean)"))
beeswarm(c(as.numeric(fitMat[metIdx,"MetFit"]), mOtherMet[metIdx]) ~ rep(1:2, rep(length(metIdx), 2)),
         add=T, cex=1.5, pwpch=pchSet[rep(fitMat[metIdx,"Met"], 2)],
         pwcol=colSet[rep(fitMat[metIdx,"Met"], 2)], corral="wrap")
dev.off()

## Primaries
mOtherPrim <- mOtherMet
pdf(paste("plots/prims_VS_otherPrims_",lungtype,"_",model,".pdf", sep=""))
primIdx <- which(!is.na(fitMat[,"PrimFit"]))
wt <- wilcox.test(as.numeric(fitMat[primIdx,"PrimFit"]), mOtherPrim[primIdx], paired=T)
boxplot(as.numeric(fitMat[primIdx,"PrimFit"]), mOtherPrim[primIdx], ylab="Normalised fitness",
        sub=sprintf("p=%.3f (paired Wilcoxon test)", wt$p.value), main="Primary cases",
        outline=F, ylim=c(0, max(c(as.numeric(fitMat[primIdx,"PrimFit"]), mOtherPrim[primIdx]))),
        names=c("Specific landscape", "Other landscapes (mean)"))
beeswarm(c(as.numeric(fitMat[primIdx,"PrimFit"]), mOtherPrim[primIdx]) ~ rep(1:2, rep(length(primIdx), 2)),
         add=T, cex=1.5, pwpch=pchSet[rep(fitMat[primIdx,"Prim"], 2)],
         pwcol=colSet[rep(fitMat[primIdx,"Prim"], 2)], corral="wrap")
dev.off()
