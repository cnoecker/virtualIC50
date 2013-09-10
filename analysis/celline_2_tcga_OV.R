source("analysis/data_utils_2.R")
source("analysis/JGLibrary.R")

library(synapseClient)
library(cgdsr)
library(glmnet)
library(parallel)
library(ROCR)
library(impute)

synapseLogin("justin.guinney@sagebase.org",'marley')

#e <- loadEntity("syn1446262")
#ov.rnaseq <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t")

e <- loadEntity("syn1687638")
ov.gistic <- as.matrix(read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t",as.is=TRUE)[,c(-1,-2)])
colnames(ov.gistic) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(ov.gistic))

e <- loadEntity("syn1446252")
ov.agilent <- as.matrix(read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t",na.strings=c("NA","null")))
is.tumor <- as.numeric(gsub("TCGA\\.\\w{2}\\.\\w{4}\\.(\\d{2}).*","\\1",colnames(ov.agilent))) < 10
is.tumor <- is.tumor & !is.na(is.tumor)
ov.agilent.g <- ov.agilent[, is.tumor]
colnames(ov.agilent.g) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(ov.agilent.g))
ov.agilent.g <- impute.knn(ov.agilent.g)$data

#is.tumor <- as.numeric(gsub("TCGA\\.\\w{2}\\.\\w{4}\\.(\\d{2}).*","\\1",colnames(ov.rnaseq))) < 10
#tmp <- ov.rnaseq[, is.tumor]
#m <- apply(tmp, 1, mean)
#tmp <- tmp[m > 1,]
#tmp <- log(tmp + 1)
#genes <- gsub("(.*?)\\|.*","\\1",rownames(tmp))
#mask <- genes != "?" & !duplicated(genes)
#ov.rnaseq.g <- tmp[mask,]
#rownames(ov.rnaseq.g) <- genes[mask]
#colnames(ov.rnaseq.g) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(tmp))

idxs <- match(colnames(ov.agilent.g), colnames(ov.gistic))
ov.agilent.m <- ov.agilent.g[, !is.na(idxs)]
ov.gistic.m <- ov.gistic[, na.omit(idxs)]

mycgds = CGDS("http://www.cbioportal.org/public-portal/")
cases <- getCaseLists(mycgds, cancerStudy="ov_tcga")
profiles <- getGeneticProfiles(mycgds,cancerStudy="luad_tcga")

muts <- getProfileData(mycgds,c("BRCA1","BRCA2","ATM","ATR","FANCA","FANCC","FANCI","FANCL","FANCD2","FANCE","FANCG","FANCM"),"ov_tcga_mutations","ov_tcga_sequenced")
muts[muts=="NaN"] = NA

e <- loadEntity("syn416789")
ov.rppa <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t")
ov.rppa <- ov.rppa[, grepl("TCGA",colnames(ov.rppa))]
colnames(ov.rppa) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(ov.rppa))

e <- loadEntity("syn1446135")
ov.clinical <- read.table(paste(e$cacheDir,e$files,sep="/"),row.names=1,header=TRUE,comment="",quote="",sep="\t",as.is=T)
ov.followup <- read.table("./resources/nationwidechildrens.org_clinical_follow_up_v1.0_ov.txt",sep="\t",header=T,quote="",as.is=T,fill=T)
uids <- gsub("(TCGA-\\w{2}-\\w{4}).*","\\1", ov.followup$bcr_followup_barcode)
ov.followup <- ov.followup[!duplicated(uids),]
ov.followup$bcr_patient_barcode <- uids[!duplicated(uids)]
ov.drug <- read.table("./resources/nationwidechildrens.org_clinical_drug_ov.txt",sep="\t",header=T,quote="",as.is=T,fill=T)

plat.sens <- getPlatinumSensitivity(ov.followup, ov.drug)
####################################

sanger <- getSanger_MetaGenomics()

common <- intersect(rownames(ov.agilent.m), featureNames(sanger))
ov.m <- ov.agilent.m[match(common, rownames(ov.agilent.m)),]
sanger.m <- sanger[match(common,featureNames(sanger)),]
tissues <- getTissueType(sampleNames(sanger.m))
ov_breast <- tissues %in% c("OVARY","BREAST")
ov_mask <- tissues %in% c("OVARY")

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
  ov.scaled <- normalize_to_X(m, sd, ov.m)
  
  y_hats <- sapply(fits, function(fit){
    predict(fit, t(ov.scaled))
  })
  y_hat <- rowMeans(y_hats)
  names(y_hat) <- colnames(ov.m)
  y_hat
}

########################################
y_hat_parp <- apply.drug.model("AZD-2281",ov_mask)
y_hat_plat_global <- apply.drug.model("Cisplatin")
y_hat_plat_ov <- apply.drug.model("Cisplatin", ov_mask)

# cisplatin analysis
idxs <- match(gsub("\\.","-",names(y_hat_plat_global)), plat.sens[,1])
y_hat_plat_global.m <- y_hat_plat_global[!is.na(idxs)]
y_hat_plat_ov.m <- y_hat_plat_ov[!is.na(idxs)]
plat.sens.m <- plat.sens[na.omit(idxs),]
surv.plat <- Surv(plat.sens.m$ttrMonths, plat.sens.m$ttrStat)

q <- quantile(y_hat_plat_global.m,c(.25, .75))
q.mask <- y_hat_plat_global.m < q[1] | y_hat_plat_global.m > q[2]
plot(survfit(surv.plat[q.mask] ~ y_hat_plat_global.m[q.mask] > median(y_hat_plat_global.m)),col=c("red","blue"),
     main="Cisplatin sensitivity model")
pval <- 1 - pchisq(foo$chisq,1)
mtext(paste("p=",format(pval,digits=2),sep=""),3)
legend("topright",fill=c("red","blue"),legend=c("Plat resistant","Plat sensitive"))


#survival analysis_
library(survival)
days <- ov.clinical$days_to_death
days[is.na(days)] <- ov.clinical$days_to_last_followup[is.na(days)]
surv <- Surv(days, !is.na(ov.clinical$days_to_death))
na.mask <- !is.na(surv)

idxs <- match(names(y_hat_parp[na.mask]), gsub("-",".",rownames(ov.clinical)[na.mask]))
coxph(surv[na.mask,][na.omit(idxs),] ~ y_hat_parp[na.mask][!is.na(idxs)])
surv.m <- surv[na.mask,][na.omit(idxs),]
yhat_parp.m <- y_hat_parp[na.mask][!is.na(idxs)]
survdiff(surv.m ~ yhat_parp.m > median(yhat_parp.m))
plot(survfit(surv.m ~ yhat_parp.m > median(yhat_parp.m)))
q <-quantile(yhat_parp.m,c(.25, .72))
quartile.mask <- yhat_parp.m < q[1] | yhat_parp.m > q[2]
f <- surv.m[quartile.mask,] ~ yhat_parp.m[quartile.mask] > median(yhat_parp.m[quartile.mask])
plot(survfit(f))
survdiff(f)

# mutation sanity checks from cbio
idxs <- match(names(y_hat_parp), rownames(muts))
y_hat_parp.m <- y_hat_parp[!is.na(idxs)]
muts.m <- muts[na.omit(idxs),]
brca.mut <- !is.na(muts.m$BRCA1) | !is.na(muts.m$BRCA2)
boxplot(y_hat_parp.m ~ brca.mut)

# mutation santiy check from tcga pub
pub_brca_muts <- read.table("../AZBrca/resources/brca_mut_tcga_paper.txt",sep="\t",header=T,as.is=T)
pub_emsy_pten <- read.table("../AZBrca/resources/emsy_pten_cna.txt",sep="\t",header=T,as.is=T)
hr.path_brcaonly <- rep(F, length(y_hat_parp))
hr.path_brca_emsy_pten <- rep(F, length(y_hat_parp))
hr.path_brcaonly[names(y_hat_parp) %in% gsub("-",".",pub_brca_muts[,1])] <- T
hr.path_brca_emsy_pten[hr.path_brcaonly | names(y_hat_parp) %in% gsub("-",".",pub_emsy_pten[,1])] <- T

boxplot(y_hat_parp ~ hr.path_brcaonly,ylab="Pred parp inhib",names=c("HR null","HR hit"),main="BRCA 1/2 mut")
mtext(paste("p=",format( wilcox.test(y_hat_parp ~ hr.path_brcaonly)$p.value,digits=2),sep=""))
boxplot(y_hat_parp ~ hr.path_brca_emsy_pten,ylab="Pred parp inhib",names=c("HR null","HR hit"),,main="BRCA 1/2 mut, PTEN, EMSY")
mtext(paste("p=",format( wilcox.test(y_hat_parp ~ hr.path_brca_emsy_pten)$p.value,digits=2),sep=""))

# fanconi genes
has.mut <- apply(muts.m, 1, function(x) any(!is.na(x)))
pat.ids <- rownames(muts.m)[has.mut]

hr.path_brca_pten_fan <- rep(F, length(y_hat_parp))
hr.path_brca_pten_fan[hr.path_brca_emsy_pten | names(y_hat_parp) %in% pat.ids] <- T

par(mfrow=c(1,3))
boxplot(y_hat_parp ~ hr.path_brcaonly,ylab="Pred parp inhib",names=c("HR null","HR hit"),main="BRCA 1/2 mut")
mtext(paste("p=",format( wilcox.test(y_hat_parp ~ hr.path_brcaonly)$p.value,digits=2),sep=""))
boxplot(y_hat_parp ~ hr.path_brca_emsy_pten,ylab="Pred parp inhib",names=c("HR null","HR hit"),,main="BRCA 1/2 mut, PTEN, EMSY")
mtext(paste("p=",format( wilcox.test(y_hat_parp ~ hr.path_brca_emsy_pten)$p.value,digits=2),sep=""))
boxplot(y_hat_parp ~ hr.path_brca_pten_fan,ylab="Pred parp inhib",names=c("HR null","HR hit"),,main="BRCA 1/2,PTEN, EMSY,FAN")
mtext(paste("p=",format( wilcox.test(y_hat_parp ~ hr.path_brca_pten_fan)$p.value,digits=2),sep=""))


#### gistic
idxs <- groupMatch(colnames(ov.m), rownames(muts),colnames(ov.gistic))
y_hat.m <- y_hat_parp[idxs[[1]]]
muts.m <- muts[idxs[[2]],]
gistic.m <- ov.gistic[,idxs[[3]]]


gistic.amp <- sort(apply(t(gistic.m), 2, function(x){
  f <- factor(x >1)
  if(length(levels(f)) < 2){ return(NA)}
  wilcox.test(y_hat.m ~ f)$p.value
}))



### rppa
idxs <- groupMatch(names(y_hat_parp), colnames(ov.rppa), rownames(muts))
y_hat.m <- y_hat_parp[idxs[[1]]]
ov.rppa.m <- ov.rppa[,idxs[[2]]]
muts.m <- muts[idxs[[3]],]

pval <- sort(p.adjust(apply(ov.rppa.m, 1, function(x){
  cor.test(x, y_hat.m, method="spearman")$p.value
}),method="BH"))



plot(as.numeric(ov.rppa.m["NF2-R-C",]), y_hat.m, ylab="olaparib pred",xlab="NF2-R-C",
     col=c("black","red")[factor(names(y_hat.m) %in% gsub("-",".",pub_brca_muts[,1]))],pch=19,cex=.9)
lines(lowess(y_hat.m ~ as.numeric(ov.rppa.m["NF2-R-C",]),f=.9),lty=2)
mtext(paste("adj pval=",format(pval["NF2-R-C"],digits=2),sep=""),3)
legend("bottomright",fill=c("black","red"),legend=c("BRCA wt","BRCA mut"))

