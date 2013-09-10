################
## colorectal models, including MSI/MSS and cimp high/low
##
source("analysis/data_utils_2.R")
source("analysis/celline_2_tcga_pipeline.R")
library(synapseClient)

synapseLogin("justin.guinney@sagebase.org",'marley')

#sanger <- getSanger_MetaGenomics()
ccle <- getCCLE_MetaGenomics()


###############
load("~/data/TCGA_ds.rda")
crc.annot <- read.table("data~/crc/TCGA_CRC_annotation.txt",sep="\t",as.is=T,header=T)
fMat <- build_feature_matrix(crc, gene.dict="cosmic",with.rppa=FALSE)
#idxs <- groupMatch(colnames(fMat), gsub("-",".",crc.annot$patient))
#fMat.m <- fMat[, idxs[[1]]]
#crc.annot.m <- crc.annot[idxs[[2]],]

#mask <- !(is.na(crc.annot.m$methylation_subtype)  | is.na(crc.annot.m$MSI_status) | is.na(crc.annot.m$hypermutated))
#fMat.m <- fMat.m[, mask]
#crc.annot.m <- crc.annot.m[mask,]
#fMat.m <- rbind(fMat.m, 
#                cimp=as.numeric(grepl("CIMP",crc.annot.m$methylation_subtype)),
#                mss=as.numeric(grepl("MSS",crc.annot.m$MSI_status)),
#                hmut=crc.annot.m$hypermutated)
                                

carcinoma_mask <- !(getTissueType(sampleNames(ccle)) %in% c("CENTRAL_NERVOUS_SYSTEM","SKIN","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
#large_intestine_mask <- getTissueType(sampleNames(ccle)) %in% "LARGE_INTESTINE"
ccle.drugs <- colnames(pData(ccle))
vIC50_2.mat <- sapply(as.list(ccle.drugs), function(drug){
  pdf(paste("results/crc/",drug,"_cosmic_perf.pdf",sep=""),width=10,height=5,useDingbats=F)
  par(mfrow=c(1,2))
  y_hat <- virtual_ic50(ccle[,carcinoma_mask],drug,crc,reverseDrug=FALSE)
  drug_fmat <- find_drug_features(y_hat,fMat, rep(TRUE, nrow(fMat)), beta_threshold=10^-3)
  dev.off()
  pdf(paste("results/crc/",drug,"_cosmic.pdf",sep=""),width=10,height=6,useDingbats=F)
  plot_features_2(drug_fmat, paste(drug," in CRC",sep=""),top=25,text.cex=.9)
  dev.off()
  y_hat
})

y_hat <- virtual_ic50(ccle[,carcinoma_mask],"RAF265",crc,reverseDrug=FALSE)
#y_hat2 <- virtual_ic50_augmentedRidge(ccle[,carcinoma_mask],"SELUMETINIB",crc,reverseDrug=FALSE)

drug_fmat <- find_drug_features(y_hat,fMat, rep(TRUE, nrow(fMat)), beta_threshold=10^-3)
top_genes <- as.character(drug_fmat$df$genes[1:20])
idxs <- groupMatch(names(y_hat), colnames(fMat))
y_hat.m <- y_hat[idxs[[1]]]
sub_fmat <- fMat[top_genes, idxs[[2]]]

M <- t(sub_fmat)

library(rpart)
fit <- rpart(y_hat.m ~ M,method="anova",
             control=rpart.control(cp=.01,xval=10,minsplit=10))
plot(fit)
text(fit,use.n=TRUE,cex=.7)
