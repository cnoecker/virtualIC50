
source("analysis/JGLibrary.R")
library(cgdsr)
library(glmnet)
library(parallel)
library(ROCR)
library(impute)


is.tumor <- function(x){
  as.numeric(gsub("TCGA\\.\\w{2}\\.\\w{4}\\.(\\d{2}).*","\\1",x)) < 10 
}

sample.ids <- function(x){
  gsub("(TCGA\\.\\w{2}\\.\\w{4}\\.\\d{2}).*","\\1",x) 
}


build.tcga.ds <- function(geneExprId, rppaId=NULL, gisticId=NULL, cbioPrefix, isRNASeq=TRUE){
  
  if(startsWith(geneExprId,"syn")){
    e <- loadEntity(geneExprId)
    geneExpr <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t")
  }else{
    header <- read.table(geneExprId,header=TRUE,row.names=1,comment="",quote="",sep="\t",nrows=1)
    geneExpr <- read.table(geneExprId,header=FALSE,row.names=1,comment="",quote="",sep="\t",skip=2)
    colnames(geneExpr) <- colnames(header)
  }
  geneExpr <- geneExpr[, grepl("^TCGA", colnames(geneExpr))]
  is.tumor <- as.numeric(gsub("TCGA\\.\\w{2}\\.\\w{4}\\.(\\d{2}).*","\\1",colnames(geneExpr))) < 10 
  tmp <- geneExpr[, is.tumor]
  if(isRNASeq){
    m <- apply(tmp, 1, mean)
    tmp <- tmp[m > 1,]
    tmp <- log(tmp + 1)
    genes <- gsub("(.*?)\\|.*","\\1",rownames(tmp))
    mask <- genes != "?" & !duplicated(genes)
    geneExpr <- tmp[mask,]
    rownames(geneExpr) <- genes[mask]
    colnames(geneExpr) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(tmp))
  }
    
  # rppa load
  if(!is.null(rppaId)){
    if(startsWith(rppaId,"syn")){
      e <- loadEntity(rppaId)
      rppa <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t")
    }else{
      rppa <- read.table(rppaId,header=TRUE,row.names=1,comment="",quote="",sep="\t")
    }
    rppa <- rppa[, grepl("TCGA",colnames(rppa))]
    colnames(rppa) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(rppa))
  }else{
    rppa <- NULL
  }
  #gistic load
  e <- loadEntity(gisticId)
  gistic <- as.matrix(read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t",as.is=TRUE)[,c(-1,-2)])
  colnames(gistic) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(gistic))

  if(startsWith(cbioPrefix,"syn")){
    muts <- buildMutationMatrix(cbioPrefix)
  }else{
    muts <- getCBIOMutationCalls(cbioPrefix)
  }
  
  return (list(geneExpr=geneExpr, rppa=rppa,gistic=gistic,mut=muts))
}

####################################

virtual_ic50 <- function(cellLineEset, drugName, tcga.dat, seed=2013, reverseDrug=FALSE){
  
  common <- intersect(rownames(tcga.dat$geneExpr), featureNames(cellLineEset))
  geneExpr.m <- tcga.dat$geneExpr[match(common, rownames(tcga.dat$geneExpr)),]
  CL.m <- cellLineEset[match(common,featureNames(cellLineEset)),]
  
  set.seed(seed)
  drug_vec <- pData(CL.m)[,drugName]
  if(reverseDrug){
    drug_vec <- -drug_vec
  }
  mask <- !is.na(drug_vec)
  X <- exprs(CL.m[,mask])
  y = drug_vec[mask]
  
  cv <- cv.glmnet(t(X),y,alpha=.1,nfolds=5)
  fits <- mclapply(1:20, function(i){
    N <- length(y)
    idxs <- sample(N,replace=TRUE)
    fit <- glmnet(t(X[,idxs]),y[idxs],alpha=.1,lambda=cv$lambda.min)
    y_hat <- rep(NA, N)
    y_hat[-idxs] <- predict(fit, t(X[,-idxs]))
    
    list(fit=fit,y_hat=y_hat)
  },mc.cores=10,mc.set.seed=TRUE)
  
  ############
  ## evaluate performance in cell lines
  
  Ymatrix <- do.call("rbind", lapply(fits, function(x) x$y_hat))
  Ymean <- colMeans(Ymatrix, na.rm=TRUE)
  Ymed <- apply(Ymatrix, 2, median, na.rm=TRUE)
  stopifnot(!any(is.na(Ymean)))
  q <- quantile(y, c(.25, .75))
  mask <- y <= q[1] | y >= q[2]
  pred <- prediction(Ymean[mask], factor(y[mask] >= q[2]))
  perf <- performance(pred, 'auc')@y.values[[1]]
  plot(performance(pred, 'tpr','fpr'),main=paste(drugName, " (n=",length(y),")",sep=""))
  text(.4, .4, labels=paste("AUC=",format(perf,digits=2),sep=""),cex=.7,pos=4)
  
  ###########
  ## apply models on new data
  m <- apply(X, 1, mean)
  sd <- apply(X, 1, sd)
  geneExpr.scaled <- normalize_to_X(m, sd, geneExpr.m)
  
  y_hats <- sapply(fits, function(x){
    predict(x$fit, t(geneExpr.scaled))
  })
  y_hat <- rowMeans(y_hats)
  names(y_hat) <- colnames(geneExpr.m)
  y_hat
}


#################################
# combine gistic and mutation and rppa and fusion and ... for lasso model
#
find_drug_features <- function(drugvec,tcga.dat, with.rppa=FALSE,beta_threshold=10^-3,num.bootstraps=50, randomize=FALSE,
                               gene.dict=c("cbio","cosmic","vogelstein"),min.count=3){
  
  #####
  type <- match.arg(gene.dict)
  switch(type){
    vogelstein = (driver.genes <- read.table("resources/vogelstein_driver_genes.txt",sep="\t",header=T,quote="",comment="",as.is=T)[,1])
    cosmic = (driver.genes <- read.table("resources/cancer_gene_census.txt",sep="\t",header=T,quote="",comment="",as.is=T)[,1])
    cbio = (driver.genes <- read.table("resources/cbio_cancer_genes.txt",sep="\t",header=T,quote="",comment="",as.is=T)[,1])
  }
  
  if(with.rppa){
    idxs <- groupMatch(names(drugvec), colnames(tcga.dat$mut), colnames(tcga.dat$gistic), colnames(tcga.dat$rppa))
  }else{
    idxs <- groupMatch(names(drugvec), colnames(tcga.dat$mut), colnames(tcga.dat$gistic))
  }
  N <- length(idxs[[1]])
  cat("Overlapping samples: ", N, "\n")
  y_hat.m <- drugvec[idxs[[1]]]
  mutM.m <- tcga.dat$mut[,idxs[[2]]]
  mutM.m <- mutM.m[rownames(mutM.m) %in% driver.genes,]
  
  gistic.m <- tcga.dat$gistic[rownames(tcga.dat$gistic) %in% rownames(mutM.m),idxs[[3]]]
  if(with.rppa){
    rppa.m <- tcga.dat$rppa[,idxs[[4]]]    
    rownames(rppa.m) <- paste("prot_", rownames(rppa.m),sep="")
  }
  
  if(randomize){
    y_hat.m <- y_hat.m[sample(length(y_hat.m))]
  }
 
  amp.m <- t(apply(gistic.m, 1, function(x) x > 1))
  del.m <- t(apply(gistic.m, 1, function(x) x < -1))
  amp.m <- amp.m[apply(amp.m, 1, sum) > min.count,]
  del.m <- del.m[apply(del.m, 1, sum) > min.count,]
  mutM.m <- mutM.m[apply(mutM.m, 1, sum) > min.count, ]
  
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
  

  if(!randomize){ pvals <- apply(A, 1, function(x){ wilcox.test(y_hat.m ~ factor(x))$p.value}) }
  if(with.rppa){
    A <- rbind(A, rppa.m)
    if(!randomize){ pvals <- c(pvals, apply(rppa.m, 1, function(x){ cor.test(x, y_hat.m)$p.value})) }
  }
  
  if(randomize){
    pvals <- rep(NA, nrow(A))
  }
  
  fits <- mclapply(1:num.bootstraps, function(i){
    N <- length(y_hat.m)
    idxs <- sample(N,replace=TRUE)
    cv.fit <- cv.glmnet(t(A)[idxs,], y_hat.m[idxs], alpha=.9,nfolds=5)
    fit <- glmnet(t(A)[idxs,], y_hat.m[idxs], lambda=cv.fit$lambda.1se,alpha=.9)
    
    y_hat <- rep(NA, N)
    y_hat[-idxs] <- predict(fit, t(A[,-idxs]))
    
    list(fit=fit,y_hat=y_hat)
  },mc.cores=min(num.bootstraps,10),mc.set.seed=TRUE)
  
  # evaluate performance
  Ymatrix <- do.call("rbind", lapply(fits, function(x) x$y_hat))
  Ymean <- colMeans(Ymatrix, na.rm=TRUE)
  Ymed <- apply(Ymatrix, 2, mean, na.rm=TRUE)
  
  na.check <- is.na(Ymean)
  if(any(na.check)){
    warning("NAs in Ymean\n")
    Ymean <- Ymean[!na.check]
    y_hat.m <- y_hat.m[!na.check]
  }
  
  q <- quantile(y_hat.m, c(.25, .75))
  mask <- y_hat.m <= q[1] | y_hat.m >= q[2]
  pred <- prediction(Ymean[mask], factor(y_hat.m[mask] >= q[2]))
  perf <- performance(pred, 'auc')@y.values[[1]]
  rho <- cor(Ymean, y_hat.m,method="spearman")
  plot(performance(pred, 'tpr','fpr'))
  text(.4, .4, labels=paste("AUC=",format(perf,digits=2),sep=""),cex=.7,pos=4)
  text(.4, .3, labels=paste("R=",format(rho,digits=2),sep=""),cex=.7,pos=4)
  
 
  R <- do.call("cbind", lapply(fits, function(x) abs(as.numeric(x$fit$beta)) > beta_threshold))
  R.plus <- rowSums(do.call("cbind", lapply(fits, function(x) as.numeric(x$fit$beta) > 0)))
  R.neg <- rowSums(do.call("cbind", lapply(fits, function(x) as.numeric(x$fit$beta) < 0)))
  pos_freq <- R.plus / (R.plus + R.neg)
  
  idxs <- order(rowSums(R),decreasing=T)
  abberationCount <- rowSums(A)[idxs]
  abberationCount[grepl("^prot",names(abberationCount))] <- NA
  
  tmp <- data.frame(genes=rownames(A)[idxs], 
                    counts=rowSums(R)[idxs],
                    posFreq=pos_freq[idxs],
                    freqCounts=rowSums(R)[idxs]/num.bootstraps,
                    noEvents=abberationCount,
                    freqEvents=abberationCount/N,
                    pvals=pvals[idxs])
  return (list(df=tmp,N=N,metric=c(rho=rho,auc=perf)))
}

buildMutationMatrix <- function(synapseId){
  e <- loadEntity(synapseId)
  tbl <- read.table(paste(e$cacheDir,"/",e$files,sep=""),header=T,as.is=T,quote="",fill=T,sep="\t")
  
  filtered.tbl <- tbl[tbl$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation"),]
  
  samples <- sort(unique(tbl$Tumor_Sample_Barcode))
  genes <- sort(unique(tbl$Hugo_Symbol))
  M <- matrix(0, nrow=length(genes),ncol=length(samples),dimnames=list(genes,samples))
  
  idxs.list <- oneToManyOperation(genes, filtered.tbl$Hugo_Symbol)
  for(i in 1:length(idxs.list)){
    barcodes <- filtered.tbl$Tumor_Sample_Barcode[idxs.list[[i]]]
    M[i, barcodes] <- 1
  }
  return (M)
}

getCBIOMutationCalls <- function(cbioPrefix, batchSize=500){
  
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  
  cancer_genes <- read.table("./resources/cbio_cancer_genes.txt",sep="\t",header=T,as.is=T,quote="")$Gene.Symbol
  cancer_genes <- cancer_genes[cancer_genes != ""]
  N <- length(cancer_genes)
  ncalls <- as.integer(N / 500) + 1
  R <- do.call("cbind", lapply(1:ncalls, function(i){
    start <- (i-1) * batchSize + 1
    end <- min(start + batchSize - 1, N)
    getProfileData(mycgds,cancer_genes[start:end],paste(cbioPrefix,"_tcga_mutations",sep=""),paste(cbioPrefix,"_tcga_sequenced",sep=""))  
  }))
  tmp <- as.matrix(data.frame(lapply(R, as.character), stringsAsFactors=FALSE))
  tmp[tmp=="NaN"] = NA
  mutM <- t(matrix(!is.na(tmp),nrow=nrow(tmp),dimnames=dimnames(R)))
  mutM
}

plot_features <- function(F,drug, disease,top=25,text.cex=.7){
  DF <- F$df[1:top,]
  N <- nrow(DF)
  sz <- DF$freqEvents * 20
  pch <- rep(21,N)
  pch[is.na(sz)] <- 22
  sz[is.na(sz)] <- 1
  col <- apply(col2rgb(rainbow(N)), 2, function(x) {rgb(x[1]/255,x[2]/255,x[3]/255,.5)})
  
  DF$pvals[DF$pvals < 10^-20] <- 10^-20

  x <- DF$freqCounts + runif(length(DF$freqCounts), 0, .1) -.05
  y <- -log10(DF$pvals)
  par(mar=c(5,5, 2,15))
  plot(x, y,pch=pch,bg=col,cex=sz,xlim=c(min(x)-.1, max(x)+.1),ylim=c(0, max(y)+1),
       ylab="-log10(pval)",xlab="Freq of selection",main=paste(drug," in ", disease," (n=",F$N,")",sep=""))
  #mtext(paste("n=",num.samples,sep=""),3)
  par(xpd=TRUE)
  
  text.col <- c("blue","black","red")[cut(DF$posFreq,breaks=c(0,.3,.7,1),include.lowest=TRUE)]
  tmp <- legend(par()$usr[2]+.03,par()$usr[4],
                legend=paste(DF$genes," (",format(DF$noEvents/F$N * 100,digits=1,trim=TRUE),"%)",sep=""),
                pch=pch,pt.bg=col,cex=text.cex,xjust=0,text.col=text.col)
  for(i in 1:N){
    x.line <- tmp$rect$left + (tmp$text$x[i] - tmp$rect$left)/2
    lines(x=c(x.line,x[i]),y=c(tmp$text$y[i],y[i]),lwd=.5,col="gray")
  }
}


