source("analysis/data_utils_2.R")
source("analysis/JGLibrary.R")

library(synapseClient)
library(cgdsr)
library(glmnet)
library(parallel)
library(ROCR)

synapseLogin("justin.guinney@sagebase.org",'marley')

e <- loadEntity("syn1446274")
read <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="")
e <- loadEntity("syn1446195")
coad <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="")

all(rownames(read) == rownames(coad))
crc <- as.matrix(cbind(read,coad))
m <- apply(crc, 1, mean)
crc <- crc[m > 1,]
crc <- log(crc + 1)
crc <- crc[,grepl("TCGA",colnames(crc))]
pat.ids <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(crc))
genes <- gsub("(.*?)\\|.*","\\1",rownames(crc))
mask <- genes != "?"
crc.g <- combine_probes_2_gene(crc[mask,],genes[mask])

## rppa
e <- loadEntity("syn1446053")
read.rppa <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="")
e <- loadEntity("syn1446043")
coad.rppa <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="")
all(rownames(coad.rppa) == rownames(read.rppa))
crc.rppa <- cbind(read.rppa, coad.rppa)
crc.rppa <- crc.rppa[, grepl("TCGA",colnames(crc.rppa))]
colnames(crc.rppa) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(crc.rppa))

####################################
drug <- "PD0325901"
ccle <- getCCLE_MetaGenomics()

common <- intersect(rownames(crc.g), featureNames(ccle))
crc.m <- crc.g[match(common, rownames(crc.g)),]
ccle.m <- ccle[match(common,featureNames(ccle)),]
tissues <- getTissueType(sampleNames(ccle.m))
non_blood <- tissues != "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"
carcinoma <- !(tissues %in% c("CENTRAL_NERVOUS_SYSTEM","MELANOMA","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
gi_mask <- tissues == "LARGE_INTESTINE"

drug.vec <- pData(ccle.m)[,drug]
mask <- !is.na(drug.vec) & carcinoma
X <- exprs(ccle.m[,mask])
y = drug.vec[mask]
tt = tissues[mask]

cv <- cv.glmnet(t(X),y,alpha=.1,nfolds=5)
fit <- glmnet(t(X),y,alpha=.1,lambda=cv$lambda.min)

fits <- mclapply(1:20, function(i){
  N <- length(y)
  idxs <- sample(N,replace=TRUE)
  weights = rep(1, N)
  weights[tt[idxs] == "LARGE_INTESTINE"] = 1
  #cv <- cv.glmnet(t(X[,idxs]),y[idxs],alpha=.1,nfolds=5)
  fit <- glmnet(t(X[,idxs]),y[idxs],weights=weights,alpha=.1,lambda=cv$lambda.1se)
  fit
  },mc.cores=10,mc.set.seed=TRUE)

m <- apply(X, 1, mean)
sd <- apply(X, 1, sd)
crc.scaled <- normalize_to_X(m, sd, crc.m)

y_hats <- sapply(fits, function(fit){
  predict(fit, t(crc.scaled))
})
y_hat <- rowMeans(y_hats)

gsets <- load.gmt.data("./resources/c2.cp.v3.1.symbols.gmt")
checkMut <- function(x){ !is.na(x) & x != "NaN" }

mycgds = CGDS("http://www.cbioportal.org/public-portal/")
cases <- getCaseLists(mycgds, cancerStudy="coadread_tcga")
profiles <- getGeneticProfiles(mycgds,cancerStudy="coadread_tcga")

gs1 <- gsets[["BIOCARTA_MAPK_PATHWAY"]]
gs2 <- gsets[["KEGG_MAPK_SIGNALING_PATHWAY"]]
muts <- getProfileData(mycgds,c("KRAS","NRAS","BRAF","PIK3R1","PIK3CA","ERBB2","ERBB3","PTEN","CTNNB1","MYC","LKB1"),"coadread_tcga_mutations","coadread_tcga_all")
mapk.muts <- getProfileData(mycgds,gs1,"coadread_tcga_mutations","coadread_tcga_all")
gistic <- getProfileData(mycgds, c("KRAS","NRAS","BRAF","ERBB2","FGFR1","IGF2","MYC","PIK3CA","CTNNB1"), "coadread_tcga_gistic","coadread_tcga_all")
mapk.gistic <- getProfileData(mycgds, gs1, "coadread_tcga_gistic","coadread_tcga_all")

muts[muts=="NaN"] = NA
mapk.muts[mapk.muts=="NaN"] <- NA
idxs <- match(pat.ids, rownames(muts))
muts.m <- muts[idxs,]
gistic.m <- gistic[idxs,]
mapk.gistic.m <- mapk.gistic[idxs,]
#rppa.m <- rppa[idxs,]
mapk.muts.m <- mapk.muts[idxs,]

has.ras.mut.nonrppa <- apply(muts.m[,c("KRAS","NRAS","BRAF")], 1, function(x){any(!is.na(x))})

common <- intersect(pat.ids, colnames(crc.rppa))
crc.rppa.m <- crc.rppa[,match(common, colnames(crc.rppa))]
y_hat.m <- y_hat[match(common, pat.ids)]
mut.m <- muts[match(common, rownames(muts)),]
has.ras.mut <- apply(mut.m[,c("KRAS","NRAS","BRAF")], 1, function(x){any(!is.na(x))})

pval <- wilcox.test(y_hat.m ~ factor(has.ras.mut))$p.value
pval.rppa <- sort(p.adjust(apply(t(crc.rppa.m), 2, function(x){
  fit <- lm(y_hat.m ~ factor(has.ras.mut) + x)
  summary(fit)$coefficients[3,4]
}),method="BH"))
                   
 pval.rppa2 <- sort(p.adjust(apply(t(crc.rppa.m), 2, function(x){
   cor.test(y_hat.m[!has.ras.mut],x[!has.ras.mut],method="spearman")$p.value
 }),method="BH"))

pval.rppa3 <- sort(p.adjust(apply(t(crc.rppa.m), 2, function(x){
  cor.test(y_hat.m,x,method="spearman")$p.value
}),method="BH"))


mask <- y_hat.m > 1
par(mfrow=c(1,2))
boxplot(y_hat ~ factor(has.ras.mut.nonrppa),names=c("RAS WT","RAS MUT"),ylab="MEK pred",outline=FALSE)
mtext(paste("pval=",format(wilcox.test(y_hat ~ factor(has.ras.mut.nonrppa))$p.value,digits=2)))
pred <- prediction(y_hat, factor(has.ras.mut.nonrppa))
perf <- performance(pred, 'auc')@y.values[[1]]
plot(performance(pred, 'tpr','fpr'))
text(.4, .4, labels=paste("AUC=",format(perf,digits=2),sep=""),cex=.7,pos=4)

plot(t(crc.rppa.m)[mask,"Fibronectin-R-C"],y_hat.m[mask],col=c("black","red")[factor(has.ras.mut)[mask]],
     pch=16,cex=.7,ylab="MEK pred",xlab="Fibronectin-R-C")
lines(lowess(t(crc.rppa.m)[mask,"Fibronectin-R-C"],y_hat.m[mask]),lty=2)
mtext(paste("adj pval=",format(pval.rppa3["Fibronectin-R-C"],digits=2),sep=""),3)
plot(t(crc.rppa.m)[mask,"AMPK_pT172-R-V"],y_hat.m[mask],col=c("black","red")[factor(has.ras.mut)[mask]],pch=16,cex=.7,ylab="MEK pred",xlab="AMPK_pT172-R-V")
lines(lowess(t(crc.rppa.m)[mask,"AMPK_pT172-R-V"],y_hat.m[mask]),lty=2)
mtext(paste("adj pval=",format(pval.rppa3["AMPK_pT172-R-V"],digits=2),sep=""),3)
plot(t(crc.rppa.m)[mask,"HER3-R-V"],y_hat.m[mask],col=c("black","red")[factor(has.ras.mut)[mask]],pch=16,cex=.7,ylab="MEK pred",xlab="HER3-R-V")
lines(lowess(t(crc.rppa.m)[mask,"HER3-R-V"],y_hat.m[mask]),lty=2)
mtext(paste("adj pval=",format(pval.rppa3["HER3-R-V"],digits=2),sep=""),3)

plot(t(crc.rppa.m)[!has.ras.mut,"beta-Catenin-R-V"],y_hat.m[!has.ras.mut],pch=16,cex=.7,ylab="MEK pred",xlab="beta-Catenin-R-V")
lines(lowess(t(crc.rppa.m)[!has.ras.mut,"beta-Catenin-R-V"],y_hat.m[!has.ras.mut]),lty=2)
mtext(paste("adj pval=",format(pval.rppa2["beta-Catenin-R-V"],digits=2),sep=""),3)

#############
## cnv/mut analysis
muts[muts=="NaN"] = NA
mapk.muts[mapk.muts=="NaN"] <- NA
idxs <- match(pat.ids, rownames(muts))
muts.m <- muts[idxs,]
gistic.m <- gistic[idxs,]
mapk.gistic.m <- mapk.gistic[idxs,]
rppa.m <- rppa[idxs,]
mapk.muts.m <- mapk.muts[idxs,]

has.ras.mut <- apply(muts.m[,c("KRAS","NRAS","BRAF")], 1, function(x){any(!is.na(x))})

mut.v <- apply(mapk.muts.m,1, function(x) any(!is.na(x)))
mut.v <- apply(mapk.gistic.m,1, function(x) any(x > 1))

pval <- wilcox.test(y_hat ~ factor(has.ras.mut))$p.value
pval2 <- wilcox.test(y_hat ~ factor(has.ras.mut | mut.v))$p.value
diff.mut <- sort(sapply(colnames(muts.m), function(gene) {
  f <- factor(!is.na(muts.m[,gene][!has.ras.mut]))
  if(length(levels(f)) == 1){ return(NA)}
  else{
    pval <- t.test(y_hat[!has.ras.mut] ~ f)$p.value
    return (pval)
  }
}))

diff.amp <- sort(apply(mapk.gistic.m, 2, function(x) {
  mask <- !is.na(x)
  pval1 <- wilcox.test(y_hat[mask] ~ factor(has.ras.mut[mask]))$p.value
  pval2 <- wilcox.test(y_hat[mask] ~ factor(has.ras.mut[mask] | x[mask] > 1))$p.value
  log(pval2) - log(pval1)
}))


wilcox.test(y_hat ~ factor(has.ras.mut | !is.na(muts.m$CTNNB1)))


pred <- prediction(y_hat, factor(has.ras.mut))
perf <- performance(pred, 'auc')@y.values[[1]]
plot(performance(pred, 'tpr','fpr'))
text(.4, .4, labels=paste("AUC=",format(perf,digits=2),sep=""),cex=.7,pos=4)