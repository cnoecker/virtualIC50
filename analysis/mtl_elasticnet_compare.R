#### #############################
###### comparison of approaches
source("analysis/data_utils_2.R")
source("analysis/celline_2_tcga_pipeline.R")
source("code/kernelMTL/bemkl_multitask_supervised_regression_variational_train.R")
source("code/kernelMTL/bemkl_multitask_supervised_regression_variational_test.R")
source("code/kernelMTL/virtual_ic50_mtl.R")

library(synapseClient)

synapseLogin("justin.guinney@sagebase.org",'marley')

sanger <- getSanger_MetaGenomics()
ccle <- getCCLE_MetaGenomics()
load("~/data/TCGA_ds.rda")


runTest <- function(disease, tcgaData, drug,iterations=10,cores=cores){
  histology <- getTissueType(sampleNames(ccle))
  match_idxs <- which(histology %in% c(disease) & !is.na(pData(ccle)[,drug]))
  N <- length(match_idxs)
  
  mtlparameters <- defaultParameters()
  
  rhos <- mclapply(1:iterations, function(i){
    test.idxs <- match_idxs[sample(N, round(N*.5))]
    
    #original algorithm
    Xtest <- exprs(ccle[,test.idxs])
    yhats <- virtual_ic50(ccle[,-test.idxs], drug,list(tcgaData[["geneExpr"]], Xtest), seed=2013,perf.eval=FALSE,num.bootstraps=3)
    #r <- virtual_ic50_mtl(ccle[,-test.idxs],drug,tcgaData,reverseDrug=FALSE, tcga.disease=disease, other.tests=list(Xtest))
    
    #### mtl 
    
    Kexpr <- rbfGeneExprKernel(list(exprs(ccle[, -test.idxs])), list(tcgaData$geneExpr), list(exprs(ccle[, test.idxs]), tcgaData$geneExpr))
    tcga.hist.factor <- rep(disease, ncol(tcgaData$geneExpr))
    Khistology <- factorKernel(list(histology[-test.idxs]), list(tcga.hist.factor), list(histology[test.idxs], tcga.hist.factor))
    
    K <- combineKernels(Kexpr, Khistology)
    vr <- virtual_ic50_mtl_v2(K$train, pData(ccle)[-test.idxs,drug],K$test, mtlparameters)
    
    # assess difference
    el_y <- yhats[[2]]
    mtl_y <- vr$y_hats[[1]]
    
    rho1 <- cor(pData(ccle)[test.idxs,drug], el_y, method="spearman")
    rho2 <- cor(pData(ccle)[test.idxs,drug], mtl_y, method="spearman")
    return (c(rho1, rho2))
  },mc.cores=cores,mc.set.seed=TRUE,mc.preschedule=FALSE)
  rhos
}

rhos_skin_raf <- runTest("SKIN",skcm,"PLX4720",iterations=10,cores=5)
boxplot(do.call("rbind",rhos_skin_raf))
rhos_lung_egfr <- runTest("LUNG",luad,"ERLOTONIB",iterations=10,cores=5)
boxplot(do.call("rbind",rhos_lung_egfr))
rhos_crc_mek <- runTest("LARGE_INTESTINE",crc,"SELUMETINIB",iterations=10,cores=5)
boxplot(do.call("rbind",rhos_crc_mek))
rhos_skin_mek <- runTest("SKIN",skcm,"SELUMETINIB",iterations=10,cores=5)
boxplot(do.call("rbind",rhos_skin_mek))

############################
########### version 2
parameters <- defaultParameters()

histology <- gsub(".*?_(.*)","\\1",rownames(pData(ccle)))
Kexpr <- rbfGeneExprKernel(list(exprs(ccle[, -test.idxs])), list(skcm$geneExpr), list(exprs(ccle[, test.idxs]), skcm$geneExpr))
tcga.hist.factor <- rep("SKIN", ncol(skcm$geneExpr))
Khistology <- factorKernel(list(histology[-test.idxs]), list(tcga.hist.factor), list(histology[test.idxs], tcga.hist.factor))
  
K <- combineKernels(Kexpr, Khistology)
#Ktrain <- abind(Kexpr[[1]], Khistology[[1]],along=3)
#Ktest <- lapply(1:length(Kexpr[[2]]), function(i) { abind(Kexpr[[2]][[i]], Kdisease[[2]][[i]],along=3)})

yhats <- virtual_ic50(ccle[,-test.idxs], "PLX4720",list(tcgaData[["geneExpr"]], exprs(ccle[,test.idxs])), seed=2013,perf.eval=FALSE,num.bootstraps=10)
vr <- virtual_ic50_mtl_v2(K$train, pData(ccle)[-test.idxs,"PLX4720"],K$test, parameters)

cor(yhats[[1]], vr$y_hats[[2]],method="spearman")
cor(yhats[[2]], vr$y_hats[[1]],method="spearman")
