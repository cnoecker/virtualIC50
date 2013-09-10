source("analysis/JGLibrary.R")
source("analysis/data_utils_2.R")
source("analysis/celline_2_tcga_pipeline.R")
source("analysis/miscFunctions.R")
source("analysis/cellline_2_tcga_mutationImporter.R")
library(synapseClient)

synapseLogin("justin.guinney@sagebase.org",'marley')

sanger <- getSanger_MetaGenomics()
ccle <- getCCLE_MetaGenomics()


###############
## build data sets
cosmicProteinPositions <- processCosmicMutationFile()

brca <- build.tcga.ds(geneExprId="syn1911247",rppaId="syn1571265",gisticId="syn1687590",cbioPrefix="brca",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
luad <- build.tcga.ds(geneExprId="syn1911220",rppaId="syn464306",gisticId="syn1687610",cbioPrefix="luad",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
coad <- build.tcga.ds(geneExprId="syn1911185",rppaId="syn1446043",gisticId="syn1687596",cbioPrefix="coadread",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
read <- build.tcga.ds(geneExprId="syn1911138",rppaId="syn1446053",gisticId="syn1687628",cbioPrefix="coadread",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
skcm <- build.tcga.ds(geneExprId="syn1910837",gisticId="syn1687618",rppaId=NULL,cbioPrefix="skcm",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
blca <- build.tcga.ds(geneExprId="syn1911172",rppaId="syn1681031",gisticId="syn1687592",cbioPrefix="blca",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
lusc <- build.tcga.ds(geneExprId="syn1917330",rppaId="syn1446049",gisticId="syn1687612",cbioPrefix="lusc",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
laml <- build.tcga.ds(geneExprId="syn1911181",rppaId=NULL,gisticId="data~/laml/laml_all_thresholded.by_genes.txt",cbioPrefix="laml",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
ucec <- build.tcga.ds(geneExprId="syn1917549",rppaId=NULL,gisticId="syn1687636",cbioPrefix="ucec",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
ov <- build.tcga.ds(geneExprId="syn1917328",rppaId=NULL,gisticId="syn1687638",cbioPrefix="ov",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
gbm <- build.tcga.ds(geneExprId="syn1911168",rppaId=NULL,gisticId="syn1687604",cbioPrefix="gbm",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
kirc <- build.tcga.ds(geneExprId="syn1911237",rppaId=NULL,gisticId="syn1687602",cbioPrefix="kirc",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)
prad <- build.tcga.ds(geneExprId="syn1917323",rppaId=NULL,gisticId="syn1687640",cbioPrefix="prad",isRNASeq=TRUE,missenseFilter=cosmicProteinPositions)



# merge rectal and colon
idxs1 <- groupMatch(rownames(coad$geneExpr), rownames(read$geneExpr))
idxs2 <- groupMatch(rownames(coad$rppa), rownames(read$rppa))
idxs3 <- groupMatch(rownames(coad$gistic), rownames(read$gistic))
crc <- list(geneExpr=cbind(coad$geneExpr[idxs1[[1]],], read$geneExpr[idxs1[[2]],]),
            rppa=cbind(coad$rppa[idxs2[[1]],], read$rppa[idxs2[[2]],]),
            gistic=cbind(coad$gistic[idxs3[[1]],], read$gistic[idxs3[[2]],]),
            mut=coad$mut)

# save datasets for fast retreival
save(brca, crc, blca, luad, skcm, lusc, laml, ucec, ov, gbm, kirc, prad, file="~/data/TCGA_ds_ver4.rda")


##########################

runDrugAnalysis <- function(drugName, drugPanelName, run.dir, drugEset, tcgaList){

  geneExprList <- lapply(tcgaList, function(x){ x$geneExpr})
  reverseDrug <- (drugPanelName == "sanger")
  
  output_prefix <- paste(run.dir,"/",drugPanelName,"_",drugName, "/",sep="")
  dir.create(output_prefix,recursive=TRUE)
  
  pdf(paste(output_prefix,"celline_perf.pdf",sep=""),width=10,height=5,useDingbats=F)
  vds <- virtual_ic50(drugEset,drug,geneExprList,reverseDrug=reverseDrug,ncores=10)
  dev.off()
  
  diseaseFMAT <- list()
  for(i in 1:length(tcgaList)){
    diseaseName <- names(tcgaList)[i]
    pdf(paste(output_prefix,diseaseName,"_perf.pdf",sep=""),width=5,height=5,useDingbats=F)
    fMat <- tcgaList[[i]]$fmat
    drug_fmat <- find_drug_features(vds$yhat[[i]]$y_mean,fMat,rep(TRUE, nrow(fMat)), beta_threshold=10^-3,ncores=10)
    dev.off()
    pdf(paste(output_prefix,diseaseName,".pdf",sep=""),width=10,height=6,useDingbats=F)
    plot_features_3(drugFMAT, paste(drug," in ",diseaseName,sep=""),top=25,text.cex=.9)
    #plot_features_2(drug_fmat, paste(drug," in ", diseaseName,sep=""),top=25,text.cex=.8)
    dev.off()
    diseaseFMAT[[diseaseName]] <- drug_fmat
  }
  vdsObj <- list(auc_celllines=vds$auc, diseaseFMAT=diseaseFMAT)
  save(vdsObj,file=paste(output_prefix,"output.rda",sep=""))
}

replotFeatures <- function(dir){
  files <- list.files(dir,)
  for(f in files){
    tmp <- strsplit(f,"_",fixed=TRUE)
    tmp1 <- tmp[[1]][1]
    drug <- tmp[[1]][2]
    cat(tmp1, " ", drug, "\n")
    load(paste(dir,f,"output.rda",sep="/"))
    diseases <- names(vdsObj$diseaseFMAT)
    for(d in diseases){
      pdf(paste(dir,"/",f,"/",d,".pdf",sep=""),width=10,height=6,useDingbats=F)
      plot_features_3(vdsObj$diseaseFMAT[[d]], paste(drug," in ",diseaseName,sep=""),top=25,text.cex=.9)
      dev.off()
    }
  }
}

replotFeatures("./results//all_2013_08_03")
################################################
# RUN ANALYSIS OVER ALL DRUGS AND DISEASES
################################################
sangerNoLymph <- sanger[,!(getTissueType(sampleNames(sanger)) %in% c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))]
ccleNoLymph <- ccle[,!(getTissueType(sampleNames(ccle)) %in% c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))]

env <- new.env()
load("~/data/TCGA_ds_ver4.rda", envir=env)
tcgaList <- as.list(env)

for(tt in names(tcgaList)){ tcgaList[[tt]]$fmat <- build_feature_matrix(tcgaList[[tt]],with.rppa=FALSE) }
lapply(tcgaList, function(x){ length(intersect(colnames(x$geneExpr), colnames(x$fmat)))})

dir.prefix <- "./results/all_2013_08_03/"

for(drug in colnames(pData(ccleNoLymph))){
  runDrugAnalysis(drug, "ccle", dir.prefix, ccleNoLymph, tcgaList)
}
for(drug in colnames(pData(sangerNoLymph))){
  runDrugAnalysis(drug, "sanger", dir.prefix, sangerNoLymph, tcgaList)
}


#################################################
# Disease specific analysis
#################################################



dir.prefix <- "results/brca/jul_12_2013/"
dir.create(dir.prefix)

prefix = "ccle"
drugEset <- ccle
reverseDrug = FALSE
tcgaDS <- brca
ncores=7

fMat <- build_feature_matrix(tcgaDS, with.rppa=FALSE)
drugs <- colnames(pData(drugEset))
#carcinoma_mask <- !(getTissueType(sampleNames(drugEset)) %in% c("CENTRAL_NERVOUS_SYSTEM","SKIN","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
carcinoma_mask <- !(getTissueType(sampleNames(drugEset)) %in% c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))

vIC50.mat <- sapply(as.list(drugs), function(drug){
  cat(drug,"\n")
  pdf(paste(dir.prefix,prefix,"_",drug,"_perf.pdf",sep=""),width=10,height=5,useDingbats=F)
  vIC50 <- virtual_ic50(drugEset[,carcinoma_mask],drug,list(tcgaDS$geneExpr),
                        reverseDrug=reverseDrug,diseaseUpweight=NULL,diseaseOffset=FALSE,ncores=ncores)
  drugFMAT <- find_drug_features(vIC50$yhat, fMat, rep(TRUE, nrow(fMat)),num.bootstraps=100,beta_threshold=10^-3,ncores=ncores)
  dev.off()
  
  pdf(paste(dir.prefix,prefix,"_",drug,".pdf",sep=""),width=10,height=6,useDingbats=F)
  plot_features_3(drugFMAT, drug,top=25,text.cex=.9)
  dev.off()
  
  save(vIC50, drugFMAT, file=paste(dir.prefix,prefix,"_",drug,".rdata",sep=""))
  
  return (NA)
})
