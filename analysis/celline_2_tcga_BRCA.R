source("analysis/data_utils_2.R")
source("analysis/JGLibrary.R")

library(synapseClient)
library(cgdsr)
library(glmnet)
library(parallel)
library(ROCR)
library(impute)

synapseLogin("justin.guinney@sagebase.org",'marley')

e <- loadEntity("syn1446183")
brca.rnaseq <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t")

is.tumor <- as.numeric(gsub("TCGA\\.\\w{2}\\.\\w{4}\\.(\\d{2}).*","\\1",colnames(brca.rnaseq))) < 10
tmp <- brca.rnaseq[, is.tumor]
m <- apply(tmp, 1, mean)
tmp <- tmp[m > 1,]
tmp <- log(tmp + 1)
genes <- gsub("(.*?)\\|.*","\\1",rownames(tmp))
mask <- genes != "?" & !duplicated(genes)
brca.rnaseq.g <- tmp[mask,]
rownames(brca.rnaseq.g) <- genes[mask]
colnames(brca.rnaseq.g) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(tmp))

e <- loadEntity("syn1571265")
brca.rppa <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t")
brca.rppa <- brca.rppa[, grepl("TCGA",colnames(brca.rppa))]
colnames(brca.rppa) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(brca.rppa))

e <- loadEntity("syn1687590")
brca.gistic <- as.matrix(read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t",as.is=TRUE)[,c(-1,-2)])
colnames(brca.gistic) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(brca.gistic))

mycgds = CGDS("http://www.cbioportal.org/public-portal/")
cases <- getCaseLists(mycgds, cancerStudy="brca_tcga")
profiles <- getGeneticProfiles(mycgds,cancerStudy="brca_tcga")
muts <- getProfileData(mycgds,c("BRCA1","BRCA2","PIK3CA"),"brca_tcga_mutations","brca_tcga_sequenced")
muts[muts=="NaN"] = NA


####################################

sanger <- getSanger_MetaGenomics()

common <- intersect(rownames(brca.rnaseq.g), featureNames(sanger))
brca.m <- brca.rnaseq.g[match(common, rownames(brca.rnaseq.g)),]
sanger.m <- sanger[match(common,featureNames(sanger)),]
tissues <- getTissueType(sampleNames(sanger.m))
ov_breast_mask <- tissues %in% c("OVARY","BREAST")
breast_mask <- tissues %in% c("BREAST")

apply.drug.model <- function(drug, sample.mask=NULL,seed=2013){
  
  set.seed(seed)
  drug.vec <- pData(sanger.m)[,drug]
  mask <- !is.na(drug.vec)
  if(!is.null(sample.mask)){
    mask <- mask & sample.mask
  }
  X <- exprs(sanger.m[,mask])
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
  brca.scaled <- normalize_to_X(m, sd, brca.m)
  
  y_hats <- sapply(fits, function(fit){
    predict(fit, t(brca.scaled))
  })
  y_hat <- rowMeans(y_hats)
  names(y_hat) <- colnames(brca.m)
  y_hat
}

########################################
y_hat_akt_1 <- apply.drug.model("A-443654",breast_mask)
y_hat_akt_2 <- apply.drug.model("AKT-inhibitor-VIII",breast_mask)
y_hat_pik3b <- apply.drug.model("AZD6482",breast_mask)
y_hat_pi3k_1 <- apply.drug.model("GDC0941",breast_mask)
y_hat_pi3k_2 <- apply.drug.model("NVP-BEZ235",breast_mask)
y_hat_mtor <- apply.drug.model("Temsirolimus",breast_mask)

l <- list("A-443654 (AKT)"=y_hat_akt_1,
          "AKT-inhibitor-VIII"=y_hat_akt_2,
          "AZD6482 (PIK3b)"=y_hat_pik3b,
          "GDC0941 (PIK3)"=y_hat_pi3k_1,
          "NVP-BEZ235 (PIK3)"=y_hat_pi3k_2,
          "Temsirolimus (MTOR)"=y_hat_mtor)

# mutation sanity checks from cbio
par(mfrow=c(2,3))
for(i in 1:length(l)){
  drug <- names(l)[i]
  y_hat <- -1 * l[[i]]
  idxs <- match(names(y_hat), rownames(muts))
  y_hat.m <- y_hat[!is.na(idxs)]
  muts.m <- muts[na.omit(idxs),]
  pik3ca.mut <- !is.na(muts.m$PIK3CA)
  f <- y_hat.m ~ pik3ca.mut
  pval <- wilcox.test(f)$p.value
  boxplot(f, ylab="sensitivity pred",names=c("PIK3CA WT", "PIK3CA MUT"),main=drug)
  mtext(paste("p=",format(pval,digits=2),sep=""))
}

### gistic
idxs <- groupMatch(colnames(brca.m), colnames(brca.gistic))
y_hat.m <- y_hat_mtor[idxs[[1]]]
gistic.m <- brca.gistic[,idxs[[2]]]

gistic.amp <- sort(apply(t(gistic.m), 2, function(x){
  f <- factor(x >1)
  if(length(levels(f)) < 2){ return(NA)}
  wilcox.test(y_hat.m ~ f)$p.value
}))

gistic.del <- sort(apply(t(gistic.m), 2, function(x){
  f <- factor(x < -1)
  if(length(levels(f)) < 2){ return(NA)}
  wilcox.test(y_hat.m ~ f)$p.value
}))

source("../pfizer_ov/src/R/chrPlot.R")
library(org.Hs.eg.db)
gistic.amp[is.na(gistic.amp)] <- 1
eg <- sapply(mget(names(gistic.amp), org.Hs.egSYMBOL2EG,ifnotfound=NA), function(x) { x[[1]] })
chr <- sapply(mget(eg[!is.na(eg)], org.Hs.egCHR,ifnotfound=NA), function(x){ x[[1]]})
loc <- sapply(mget(eg[!is.na(eg)], org.Hs.egCHRLOC,ifnotfound=NA), function(x){ x[[1]]})
mask <- !is.na(chr) & !is.na(loc) & chr %in% c(1:22)
foo = cbind(loc=abs(loc[mask]),chr=as.numeric(chr[mask]),gene=names(gistic.amp[!is.na(eg)][mask]),val=-log10(gistic.amp[!is.na(eg)][mask]))
amp.count <- apply(gistic.m,1,function(x){ sum(x > 1)}) / ncol(gistic.m)
pdf("foo.pdf",width=10,height=6)
par(mfrow=c(2,1))
plotGWASPValues.default(abs(loc[mask]), as.numeric(chr[mask]), -log10(gistic.amp[!is.na(eg)][mask]),cex=.5,yLabel="-log10(pval)",mainTitle="")
plotGWASPValues.default(abs(loc[mask]), as.numeric(chr[mask]), amp.count[!is.na(eg)][mask],cex=.5,yLabel="amp freq",mainTitle="")
dev.off()

gistic.del[is.na(gistic.del)] <- 1
eg <- sapply(mget(names(gistic.del), org.Hs.egSYMBOL2EG,ifnotfound=NA), function(x) { x[[1]] })
chr <- sapply(mget(eg[!is.na(eg)], org.Hs.egCHR,ifnotfound=NA), function(x){ x[[1]]})
loc <- sapply(mget(eg[!is.na(eg)], org.Hs.egCHRLOC,ifnotfound=NA), function(x){ x[[1]]})
mask <- !is.na(chr) & !is.na(loc) & chr %in% c(1:22)
foo = cbind(loc=abs(loc[mask]),chr=as.numeric(chr[mask]),gene=names(gistic.del[!is.na(eg)][mask]),val=-log10(gistic.del[!is.na(eg)][mask]))


library(ggbio)
data(hg19IdeogramCyto, package = "biovizBase")
data(hg19Ideogram, package = "biovizBase")
library(GenomicRanges)

chrs <- as.character(levels(seqnames(hg19IdeogramCyto)))
seqlths <- seqlengths(hg19Ideogram)[chrs]


gr.cna <- GRanges(chr[mask],IRanges(start=abs(loc[mask]),width=1),pval=-log10(gistic.del[!is.na(eg)][mask]))
gr.cna <- keepSeqlevels(gr.cna, as.character(1:22))
#names(seqlths) <- gsub("chr", "", names(seqlths))
#seqlengths(gr.cna) <- seqlths[names(seqlengths(gr.cna))]

plotGrandLinear(gr.cna, aes(y = pval))



### rppa
y_hat <- y_hat_mtor
idxs <- groupMatch(names(y_hat), colnames(brca.rppa))
y_hat.m <- y_hat[idxs[[1]]]
brca.rppa.m <- brca.rppa[,idxs[[2]]]

pval <- sort(p.adjust(apply(brca.rppa.m, 1, function(x){
  cor.test(x, y_hat.m, method="spearman")$p.value
}),method="BH"))


pdf("foo.pdf")
plot(as.numeric(brca.rppa.m["4E-BP1_pS65-R-V",]), y_hat.m, ylab="mtor pred",xlab="4E-BP1_pS65-R-V",,pch=19,cex=.9)
lines(lowess(y_hat.m ~ as.numeric(brca.rppa.m["4E-BP1_pS65-R-V",]),f=.9),lty=2)
mtext(paste("adj pval=",format(pval["4E-BP1_pS65-R-V"],digits=2),sep=""),3)
dev.off()


###########################
## lasso model

# mutation data from cbio
cancer_genes <- read.table("./resources/cbio_cancer_genes.txt",sep="\t",header=T,as.is=T,quote="")
cancer.muts1 <- getProfileData(mycgds,cancer_genes$Gene.Symbol[1:1000],"brca_tcga_mutations","brca_tcga_sequenced")
cancer.muts2 <- getProfileData(mycgds,cancer_genes$Gene.Symbol[1001:2000],"brca_tcga_mutations","brca_tcga_sequenced")
cancer.muts3 <- getProfileData(mycgds,cancer_genes$Gene.Symbol[2001:nrow(cancer_genes)],"brca_tcga_mutations","brca_tcga_sequenced")
all(rownames(cancer.muts1) == rownames(cancer.muts2))
all(rownames(cancer.muts1) == rownames(cancer.muts3))
cancer.muts <- cbind(cancer.muts1,cancer.muts2, cancer.muts3)
rownames(cancer.muts) <- rownames(cancer.muts1)
tmp <- as.matrix(data.frame(lapply(cancer.muts, as.character), stringsAsFactors=FALSE))
tmp[tmp=="NaN"] = NA
mutM <- t(matrix(!is.na(tmp),nrow=nrow(tmp),dimnames=dimnames(cancer.muts)))

# combine gistic and mutation data
find_drug_features <- function(drugvec,with.rppa=FALSE){
  if(with.rppa){
    idxs <- groupMatch(names(drugvec), colnames(mutM), colnames(brca.gistic), colnames(brca.rppa))
  }else{
    idxs <- groupMatch(names(drugvec), colnames(mutM), colnames(brca.gistic))
  }
  y_hat.m <- drugvec[idxs[[1]]]
  mutM.m <- mutM[,idxs[[2]]]
  brca.gistic.m <- brca.gistic[rownames(brca.gistic) %in% rownames(mutM),idxs[[3]]]
  if(with.rppa){
    brca.rppa.m <- brca.rppa[,idxs[[4]]]    
  }
  
  amp.m <- t(apply(brca.gistic.m, 1, function(x) x > 1))
  del.m <- t(apply(brca.gistic.m, 1, function(x) x < -1))
  amp.m <- amp.m[apply(amp.m, 1, sum) > 0,]
  del.m <- del.m[apply(del.m, 1, sum) > 0,]
  mutM.m <- mutM.m[apply(mutM.m, 1, sum) > 0, ]
  
  idxs <- groupMatch(rownames(mutM.m), rownames(amp.m))
  mut_amp <- mutM.m[idxs[[1]],] | amp.m[idxs[[2]],]
  idxs <- groupMatch(rownames(mutM.m), rownames(del.m))
  mut_del <- mutM.m[idxs[[1]],] | del.m[idxs[[2]],]
  
  rownames(amp.m) <- paste("amp_", rownames(amp.m),sep="")
  rownames(del.m) <- paste("del_", rownames(del.m),sep="")
  rownames(mutM.m) <- paste("mut_", rownames(mutM.m),sep="")
  rownames(mut_amp) <- paste("mut_amp_", rownames(mut_amp),sep="")
  rownames(mut_del) <- paste("mut_del_", rownames(mut_del),sep="")
  
  A <- rbind(amp.m, del.m, mutM.m, mut_amp, mut_del)
  A <- A[rowSums(A) > 2,]
  if(with.rppa){
    A <- rbind(A, brca.rppa.m)
  }
  
  fits <- mclapply(1:50, function(i){
    N <- length(y_hat.m)
    idxs <- sample(N,replace=TRUE)
    cv.fit <- cv.glmnet(t(A)[idxs,], y_hat.m[idxs], alpha=.9,nfolds=5)
    fit <- glmnet(t(A)[idxs,], y_hat.m[idxs], lambda=cv.fit$lambda.1se,alpha=.9)
    fit
  },mc.cores=10,mc.set.seed=TRUE)
  
  R <- do.call("cbind", lapply(fits, function(fit) abs(as.numeric(fit$beta)) > 10^-3))
  idxs <- order(rowSums(R),decreasing=T)
  tmp <- data.frame(genes=rownames(A)[idxs], counts=rowSums(R)[idxs],noAberrations=rowSums(A)[idxs])
  return (tmp)
}

R.mtor <- find_drug_features(y_hat_mtor)
R.azd6482 <- find_drug_features(y_hat_pik3b)
R.akt1 <- find_drug_features(y_hat_akt_1)

R.mtor.rppa <- find_drug_features(y_hat_mtor,with.rppa=T)
R.azd6482.rppa <- find_drug_features(y_hat_pik3b,with.rppa=T)
R.akt1.rppa <- find_drug_features(y_hat_akt_1,with.rppa=T)
