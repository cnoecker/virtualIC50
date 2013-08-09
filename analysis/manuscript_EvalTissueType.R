source("analysis/JGLibrary.R")
source("analysis/data_utils_2.R")
source("analysis/celline_2_tcga_pipeline.R")
source("analysis/miscFunctions.R")
source("analysis/cellline_2_tcga_mutationImporter.R")
library(synapseClient)

ncores=5

if(!exists("ccle",.GlobalEnv)){
  synapseLogin("justin.guinney@sagebase.org",'marley')
  ccle <- getCCLE_MetaGenomics()
}

if(!exists("tcgaList",.GlobalEnv)){
  env <- new.env()
  load("~/data/TCGA_ds_ver4.rda", envir=env)
  tcgaList <- as.list(env)
}

drugs <- c("ERLOTONIB","SELUMETINIB","LAPATINIB","VANDETANIB","SORAFENIB","TKI258","PLX4720")
diseases <- c("luad","lusc","skcm","ov","crc","brca","ucec")
cldiseases <- c("ALL","LUAD","LUSC","SKIN","OVARY","LARGE_INTESTINE","BREAST","ENDOMETRIUM")


virtual_ic50_SS <- function(cellLineEset, drugName, unlabeledExpr, seed=2013, reverseDrug=FALSE,perf.eval=TRUE,
                            num.bootstraps=20,ncores=5,alpha=.1){

  set.seed(seed)
  
  drug_vec <- pData(cellLineEset)[,drugName]
  if(reverseDrug){
    drug_vec <- -drug_vec
  }
  has_drug_mask <- !is.na(drug_vec)
  
  isSS <- !is.null(unlabeledExpr)
  lX <- scale(t(exprs(cellLineEset[,has_drug_mask])))
  y = drug_vec[has_drug_mask]
  tissueTypes <- getTissueType(cellLineEset)[has_drug_mask]
  N <- length(y)
  if(isSS){
    Nu <- nrow(unlabeledExpr)
    y_c <- c(y - mean(y), rep(0, Nu))
    X <- rbind(lX,unlabeledExpr)
  }else{
    y_c <- y - mean(y)
    X <- lX
  }
 
  cv <- cv.glmnet(X,y_c,alpha=alpha,nfolds=5,standardize=FALSE)
  
  fits <- mclapply(1:num.bootstraps, function(i){
    idxs <- sample(N,replace=TRUE)
    if(isSS){
      y_c <- c((y - mean(y))[idxs], rep(0, Nu))
      X <- rbind(lX[idxs,],unlabeledExpr)
    }else{
      y_c <- (y - mean(y))[idxs]
      X <- lX[idxs,]
    }
    
    fit <- glmnet(X,y_c,alpha=alpha,lambda=cv$lambda.min,standardize=FALSE)
    y_hat <- rep(NA, N)
    y_hat[-idxs] <- predict(fit, lX[-idxs,])
    list(y_hat=y_hat)
  },mc.cores=ncores,mc.set.seed=TRUE,mc.preschedule=FALSE)
  
  perfList <- list()
    
  
  Ymatrix <- do.call("rbind", lapply(fits, function(x) x$y_hat))
  Ymean <- colMeans(Ymatrix, na.rm=TRUE)
  na.mask <- !is.na(Ymean)
  Ymean <- Ymean[na.mask]
  y <- y[na.mask]
  tt <- tissueTypes[na.mask]
  ttTbl <- table(tt)
  diseaseTests <- c("ALL",names(ttTbl)[ttTbl >= 10])
  
  for(disease in diseaseTests){
    if(disease == "ALL"){
      diseaseMask <- rep(TRUE, length(y))
    }else{
      diseaseMask <- tt %in% disease
    }
    
    cat(disease,"\n")
    yDisease <- y[diseaseMask]
    yDiseaseMean <- Ymean[diseaseMask]
    q <- quantile(yDisease, c(.40, .60))
    mask <- yDisease <= q[1] | yDisease >= q[2]
    #browser()
    pred <- prediction(yDiseaseMean[mask], factor(yDisease[mask] >= q[2]))
    auc <- performance(pred, 'auc')@y.values[[1]]
    quant_80 <- yDisease > quantile(yDisease, .80)
    ppv <- sum((yDiseaseMean > quantile(yDiseaseMean, .80)) & quant_80) / sum(quant_80)
    quant_20 <- yDisease < quantile(yDisease, .20)
    npv <- sum((yDiseaseMean < quantile(yDiseaseMean, .20)) & quant_20) / sum(quant_20)
    
    rho <- cor.test(yDisease, yDiseaseMean, method="spearman")
    perfList[[disease]] <- list(auc=auc,rho=rho,npv=npv,ppv=ppv)
  }
  perfList
}


ccleNoLymph <- ccle[,!(getTissueType(sampleNames(ccle)) %in% c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))]
load("/home/cferte/resources/cell_line_infos.rda")
ccleinfo <- cells.info$ccle.info
tt <- getTissueType(sampleNames(ccleNoLymph))
histology <- ccleinfo[sampleNames(ccleNoLymph),"Hist.Subtype1"]
clnames <- gsub("(.*?)_.*","\\1", sampleNames(ccleNoLymph), perl=TRUE)
luad.mask <- tt == "LUNG" & histology == "adenocarcinoma"
lusc.mask <- tt == "LUNG" & histology == "squamous_cell_carcinoma"

lbls <- sampleNames(ccleNoLymph)
lbls[luad.mask] <- paste(clnames[luad.mask],"_LUAD",sep="")
lbls[lusc.mask] <- paste(clnames[lusc.mask],"_LUSC",sep="")
sampleNames(ccleNoLymph) <- lbls

tcgaList <- tcgaList[diseases]
common <- featureNames(ccleNoLymph)

for(i in seq_along(tcgaList)){ common <- intersect(common, rownames(tcgaList[[i]]$geneExpr))}
ccleNoLymph.m <- ccleNoLymph[match(common,featureNames(ccleNoLymph)),]
geneExprList <- lapply(tcgaList, function(tmp){ 
  x <- tmp$geneExpr
  scale(t(x[match(common,rownames(x)),]))
})

R <- list()
for(drug in drugs){
  for(d in diseases){
    Xu <- geneExprList[[d]]
    testName <- paste(d,"_",drug,sep="")
    cat(testName,"\n")
    R[[testName]] <- virtual_ic50_SS(ccleNoLymph.m, drug, Xu, seed=2013, reverseDrug=FALSE,perf.eval=TRUE,num.bootstraps=20,ncores=ncores, alpha=.1)
  }
  testName <- paste("NA_",drug,sep="")
  R[[testName]] <- virtual_ic50_SS(ccleNoLymph.m, drug, NULL, seed=2013, reverseDrug=FALSE,perf.eval=TRUE,num.bootstraps=20,ncores=ncores, alpha=.1)
}

save(R,file="./results/tissue_specific_evaluation.rda")

#tissueTypes <- c("ALL",unique(getTissueType(colnames(ccle))))
aucM <- matrix(NA, nrow=length(R),ncol=length(cldiseases),dimnames=list(names(R), cldiseases))
rhoM <- matrix(NA, nrow=length(R),ncol=length(cldiseases),dimnames=list(names(R), cldiseases))
for(i in seq_along(R)){
  for(tt in cldiseases){
    aucM[names(R)[i],tt] <- R[[i]][[tt]]$auc
  }
}
drug_aucM <- lapply(drugs, function(x) aucM[grepl(x, rownames(aucM)),])
names(drug_aucM) <- drugs

drugPlot <- function(M,drug){
  par(mfrow=c(1,7),oma = c(0, 0, 3, 0))
  for(i in 2:(ncol(M))){
    col=rep("gray",8)
    col[i-1] <- "red"
    col[8] <- "blue"
    barplot(M[,i], ylim=c(0,1),names.arg="",col=col,main=colnames(M)[i])
    abline(h=.5, lty=2)
  }
  mtext(drug,cex=1.5,outer=TRUE)
}
for(i in seq_along(drug_aucM)){
  M <- drug_aucM[[i]]
  drug <- strsplit(rownames(M)[1],"_")[[1]][2]
  pdf(paste("./results/ss_tissue_eval/",drug,".pdf",sep=""),width=10,height=4)
  drugPlot(M,drug)
  dev.off()
}
