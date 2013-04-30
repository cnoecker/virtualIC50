source("analysis/data_utils_2.R")
source("analysis/JGLibrary.R")

library(synapseClient)
library(cgdsr)
library(glmnet)
library(parallel)
library(ROCR)

synapseLogin("justin.guinney@sagebase.org",'marley')

e <- loadEntity("syn1687610")
luad.gistic <- as.matrix(read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t",as.is=TRUE)[,c(-1,-2)])
colnames(luad.gistic) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(luad.gistic))

e <- loadEntity("syn418003")
luad.rnaseq <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t")

e <- loadEntity("syn1571581")
luad.mut <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,comment="",quote="",sep="\t",as.is=T,fill=T)

is.tumor <- as.numeric(gsub("TCGA\\.\\w{2}\\.\\w{4}\\.(\\d{2}).*","\\1",colnames(luad.rnaseq))) < 10
tmp <- luad.rnaseq[, is.tumor]
m <- apply(tmp, 1, mean)
tmp <- tmp[m > 1,]
tmp <- log(tmp + 1)
genes <- gsub("(.*?)\\|.*","\\1",rownames(tmp))
mask <- genes != "?" & !duplicated(genes)
luad.rnaseq.g <- tmp[mask,]
rownames(luad.rnaseq.g) <- genes[mask]
colnames(luad.rnaseq.g) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(tmp))

idxs <- match(colnames(luad.rnaseq.g), colnames(luad.gistic))
luad.rnaseq.m <- luad.rnaseq.g[, !is.na(idxs)]
luad.gistic.m <- luad.gistic[, na.omit(idxs)]

mycgds = CGDS("http://www.cbioportal.org/public-portal/")
cases <- getCaseLists(mycgds, cancerStudy="luad_tcga")
profiles <- getGeneticProfiles(mycgds,cancerStudy="luad_tcga")

muts <- getProfileData(mycgds,c("KRAS","NRAS","HRAS","BRAF","STK11","NF1","EGFR","PIK3CA","TP53"),"luad_tcga_mutations","luad_tcga_all")
muts[muts=="NaN"] = NA

e <- loadEntity("syn464306")
luad.rppa <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t")
luad.rppa <- luad.rppa[, grepl("TCGA",colnames(luad.rppa))]
colnames(luad.rppa) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(luad.rppa))


####################################

ccle <- getCCLE_MetaGenomics()

common <- intersect(rownames(luad.rnaseq.m), featureNames(ccle))
luad.m <- luad.rnaseq.m[match(common, rownames(luad.rnaseq.m)),]
ccle.m <- ccle[match(common,featureNames(ccle)),]
tissues <- getTissueType(sampleNames(ccle.m))
carcinoma <- !(tissues %in% c("CENTRAL_NERVOUS_SYSTEM","MELANOMA","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
lung.mask <- tissues %in% "LUNG"

apply.drug.model <- function(drug, sample.mask=NULL,seed=2013){

  set.seed(seed)
  drug.vec <- pData(ccle.m)[,drug]
  mask <- !is.na(drug.vec)
  if(!is.null(sample.mask)){
    mask <- mask & sample.mask
  }
  X <- exprs(ccle.m[,mask])
  y = drug.vec[mask]
  
  cv <- cv.glmnet(t(X),y,alpha=.1,nfolds=5)
  #fit <- glmnet(t(X),y,alpha=.1,lambda=cv$lambda.min)
  
  fits <- mclapply(1:20, function(i){
    N <- length(y)
    idxs <- sample(N,replace=TRUE)
    #cv <- cv.glmnet(t(X[,idxs]),y[idxs],alpha=.1,nfolds=5)
    fit <- glmnet(t(X[,idxs]),y[idxs],alpha=.1,lambda=cv$lambda.1se)
    fit
  },mc.cores=10,mc.set.seed=TRUE)
  
  m <- apply(X, 1, mean)
  sd <- apply(X, 1, sd)
  luad.scaled <- normalize_to_X(m, sd, luad.m)
  
  y_hats <- sapply(fits, function(fit){
    predict(fit, t(luad.scaled))
  })
  y_hat <- rowMeans(y_hats)
  names(y_hat) <- colnames(luad.m)
  y_hat
}

########################################
y_hat_mek <- apply.drug.model("PD0325901",carcinoma)
y_hat_erl <- apply.drug.model("ERLOTONIB",lung.mask)

# mutation sanity checks
idxs <- match(names(y_hat_mek), rownames(muts))
y_hat_mek.m <- y_hat_mek[!is.na(idxs)]
y_hat_erl.m <- y_hat_erl[!is.na(idxs)]

muts.m <- muts[na.omit(idxs),]
egfr.mut <- !is.na(muts.m$EGFR)
ras.mut <- !is.na(muts.m$KRAS) | !is.na(muts.m$NRAS) | (muts.m$BRAF=="V600E" & !is.na(muts.m$BRAF))
stk11 <- !is.na(muts.m$STK11)
                                                                              
boxplot(y_hat_erl.m ~ !is.na(muts.m$EGFR), outline=F,main="Erlotinib",names=c("EGFR WT","EGFR MUT"),ylab="Erlotinib pred")
pval <- wilcox.test(y_hat_erl.m ~ !is.na(muts.m$EGFR))$p.value
mtext(paste("p=",format(pval,digits=2),sep=""))
                                                                              
pval <- wilcox.test(y_hat_mek.m ~ ras.mut)$p.value
boxplot(y_hat_mek.m ~ ras.mut,main="MEK",outline=F,main="MEK",names=c("RAS WT","RAS MUT"),ylab="MEK pred")
mtext(paste("p=",format(pval,digits=2)),3)


#### gistic
idxs <- groupMatch(colnames(luad.m), rownames(muts),colnames(luad.gistic))
y_hat.m <- y_hat[idxs[[1]]]
muts.m <- muts[idxs[[2]],]
gistic.m <- luad.gistic[,idxs[[3]]]

ras.mut <- !is.na(muts.m$KRAS) | !is.na(muts.m$NRAS)
wilcox.test(y_hat.m ~ ras.mut)

gistic.pvals.2 <- sort(apply(gistic.m, 1, function(x){
  f <- factor(x > 1)
  if(length(levels(f)) < 2){ return(NA)}
  wilcox.test(y_hat.m ~ f)$p.value
}))



### rppa
idxs <- groupMatch(names(y_hat_mek), colnames(luad.rppa), rownames(muts))
y_hat_mek.m <- y_hat_mek[idxs[[1]]]
y_hat_erl.m <- y_hat_erl[idxs[[1]]]
luad.rppa.m <- luad.rppa[,idxs[[2]]]
muts.m <- muts[idxs[[3]],]

pval.mek <- sort(apply(luad.rppa.m, 1, function(x){
  cor.test(x, y_hat_mek.m, method="spearman")$p.value
}))
pval.erl <- sort(apply(luad.rppa.m, 1, function(x){
  cor.test(x, y_hat_erl.m, method="spearman")$p.value
}))

plot(as.numeric(luad.rppa.m["EGFR_pY1068-R-V",]), y_hat_erl.m, ylab="ERLOTINIB pred",xlab="EGFR_pY1068-R-V",
     col=c("black","red")[factor(!is.na(muts.m$EGFR))],pch=19,cex=.9)
lines(lowess(y_hat_erl.m ~ as.numeric(luad.rppa.m["EGFR_pY1068-R-V",]),f=.9),lty=2)
mtext(paste("p=",format(pval.erl["EGFR_pY1068-R-V"],digits=2),sep=""),3)
legend("bottomright",fill=c("black","red"),legend=c("EGFR wt","EGFR mut"))

