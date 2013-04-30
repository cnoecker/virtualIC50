
source("analysis/JGLibrary.R")
library(cgdsr)
library(glmnet)
library(parallel)
library(ROCR)
library(impute)

create.methylation.library <- function(geneExprId, methylationId){
  library(IlluminaHumanMethylation450k.db)
  if(startsWith(geneExprId,"syn")){
    e <- loadEntity(geneExprId)
    geneExpr <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t")
  }else{
    header <- read.table(geneExprId,header=TRUE,row.names=1,comment="",quote="",sep="\t",nrows=1)
    geneExpr <- read.table(geneExprId,header=FALSE,row.names=1,comment="",quote="",sep="\t",skip=2)
    colnames(geneExpr) <- colnames(header)
  }
  geneExpr <- data.matrix(geneExpr[, grepl("^TCGA", colnames(geneExpr))])
  tmp <- geneExpr
  if(isRNASeq){
    m <- apply(tmp, 1, mean)
    tmp <- tmp[m > 1,]
    tmp <- log(tmp + 1)
    genes <- gsub("(.*?)\\|.*","\\1",rownames(tmp))
    mask <- genes != "?" & !duplicated(genes)
    geneExpr <- tmp[mask,]
    rownames(geneExpr) <- genes[mask]
  }
  
  
  e <- loadEntity("syn415444") #luad
 # e <- loadEntity("syn1445956") #brca
  ncol <- ncol(read.table(paste(e$cacheDir,"/",e$files,sep=""),comment.char="",header=T,row.names=1,as.is=T,,nrows=2))
  
  cols <- c("character",rep("numeric",ncol-1))
  tbl <- data.matrix(read.table(paste(e$cacheDir,"/",e$files,sep=""),comment.char="",header=T,row.names=1,as.is=T,colClasses=cols))
  # remove large numbers of NAs
  na.count <- apply(tbl, 1, function(x) sum(is.na(x)))
  threshold <- floor(ncol * .05)
  tbl <- tbl[na.count < threshold,]
  imputed.tbl <- impute.knn(tbl)$data
  
  genes <- unlist(mget(rownames(imputed.tbl), IlluminaHumanMethylation450kSYMBOL,ifnotfound=NA))
  idxs <- oneToManyOperation(rownames(geneExpr), genes)
  
  sample.idxs <- groupMatch(sample.ids(colnames(geneExpr)), sample.ids(colnames(imputed.tbl)))
  gene.expr.m <- geneExpr[, sample.idxs[[1]]]
  colnames(gene.expr.m) <- sample.ids(colnames(gene.expr.m))
  methyl.m <- imputed.tbl[, sample.idxs[[2]]]
  colnames(methyl.m) <- sample.ids(colnames(methyl.m))
  tumor.mask <- is.tumor(colnames(gene.expr.m))
  
  cat("Num tumor: ", sum(tumor.mask),"\n")
  cat("Num normal: ", sum(!tumor.mask),"\n")
  
  gene.expr.t <- gene.expr.m[rownames(gene.expr.m) %in% cancer_genes, tumor.mask]
  methyl.t <- methyl.m[, tumor.mask]
  
  idxs <- oneToManyOperation(rownames(gene.expr.t), genes)
  
  sapply(1:length(gene.expr.t), function(i){
    rhos <- cor(t(gene.expr.t[i,,drop=FALSE]), t(methyl.t[idxs[[i]],,drop=FALSE]),method="spearman"))
    which <- rhos < -.25
    minratio <- rep(NA, length(rhos))
    for(j in 1:length(rhos)){
      if(rhos[j] < -.25){
        sapply(seq(.1,.9,by=.1), function(cutoff){ 
          mask <- methyl.t[idxs[[i]][j]] < cutoff
          rho1 <- gene.expr.t[i]
        })
      }
    }
  
  })
  
  m.pipeline <- mclapply(as.list(1:nrow(gene.expr.m)), function(i){
    if(is.null(idxs[[i]])) { return (NA)}
    tcgaEpigeneticSilencingTests(gene.expr.m[i,], methyl.m[idxs[[i]],,drop=FALSE], tumor.mask)  
  },mc.cores=10)
  names(m.pipeline) <- rownames(gene.expr.m)
  
  foo <- sapply(as.list(1:nrow(gene.expr.m)), function(i){
    if(is.null(idxs[[i]])) { return (NA)}
    min(cor(t(gene.expr.m[i,tumor.mask,drop=FALSE]), t(methyl.m[idxs[[i]],tumor.mask,drop=FALSE]),method="spearman"))
  })
  m <- m.pipeline[!is.na(m.pipeline)]
  tmp <- sapply(m, function(x){ any(x$stringent > 3) } )
  plot(gene.expr.m)
}

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

  muts <- getCBIOMutationCalls(cbioPrefix)
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
find_drug_features <- function(drugvec,tcga.dat, with.rppa=FALSE,beta_threshold=10^-3,num.bootstraps=50, randomize=FALSE,gene.dict="cbio",min.count=3){
  
  #####
  stopifnot(gene.dict %in% c("cbio","cosmic","vogelstein"))
  if(gene.dict == "vogelstein"){
    driver.genes <- read.table("resources/driver_genes.txt",sep="\t",header=T,quote="",comment="",as.is=T)[,1]
  }else if(gene.dict == "cbio"){
    driver.genes <- read.table("resources/cancer_gene_census.txt",sep="\t",header=T,quote="",comment="",as.is=T)[,1]
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
  if(driver.genes.only){
    mutM.m <- mutM.m[rownames(mutM.m) %in% driver.genes,]
  }
  
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

# ##################################
# compute the 4 metrics specified in tcga ovarian paper (doi:10.1038/nature10166)
#   - geneExpr - a named vector of gene level expression with names equal to patient barcode (matching column names of cpgMethMat)
#   - cpgMethMat - a matrix of methylation beta values which
#     - has rows per cpg that correspond to gene expression measured
#     - has columns in the same order as samples in the expression vector
#   - tumor - a 0/1 vector specifying if a sample measuring tumor tissue or normal (fallopian tube)
#   - pctHyper - % samples tested for hypermethylation - default as in paper of 0.10
#   - relaxedThreshold - vector of 'relaxed' cutoffs for tests 1-4 - defaults as in paper
#   - stringentThreshold - vector of 'stringent' cutoffs for tests 1-4 - defaults as in paper
#
tcgaEpigeneticSilencingTests <- function(geneExpr, cpgMethMat, tumor, pctHyper=0.10, relaxedThreshold = c(0.5, 0.1, 2, -0.2), stringentThreshold = c(0.4, 0.3, 3, -0.3)){
  
  res <- list()
  
  ## TEST 1
  ## MEAN BETA IN NORMALS
  res$test1 <- rowMeans(cpgMethMat[, tumor==0,drop=FALSE])
  names(res$test1) <- rownames(cpgMethMat)
  
  ## TEST 2
  ## DIFFERENCE IN BETA BETWEEN MEAN NORMALS AND 90%ILE OF TUMORS
  res$test2 <- sapply(as.list(rownames(cpgMethMat)), function(x){
    quantile(cpgMethMat[x, tumor==1], probs=(1-pctHyper), na.rm=T) - mean(cpgMethMat[x, tumor==0])
  })
  names(res$test2) <- rownames(cpgMethMat)
  
  ## TEST 3
  ## EXPRESSION FC BETWEEN NORMALS AND TOP 10% TUMORS W/ HIGHEST METHYLATION (BETA)
  res$test3 <- sapply(as.list(rownames(cpgMethMat)), function(x){
    tmp <- cpgMethMat[x, tumor == 1]
    theseTumors <- names(tmp)[ tmp >= quantile(tmp, probs=(1-pctHyper), na.rm=T) ]
    expr <- c(geneExpr[theseTumors], geneExpr[tumor==0])
    tum <- c(rep(1, length(theseTumors)), rep(0, sum(tumor==0)))
    tmpFit <- lm(expr ~ tum)
    2^(-1*tmpFit$coefficients["tum"])
  })
  names(res$test3) <- rownames(cpgMethMat)
  
  ## TEST 4
  ## CORRELATION BETWEEN METHYLATION AND EXPRESSION VALUES
  res$test4 <- sapply(as.list(1:nrow(cpgMethMat)), function(x){ cor(cpgMethMat[x, ], geneExpr, method="spearman", use="complete.obs")})
  names(res$test4) <- rownames(cpgMethMat)
  
  ## RELAXED AND STRINGENT THRESHOLDS FROM TCGA PAPER
  res$relaxed <- (res$test1 < relaxedThreshold[1]) + (res$test2 > relaxedThreshold[2]) + (res$test3 > relaxedThreshold[3]) + (res$test4 < relaxedThreshold[4])
  res$stringent <- (res$test1 < stringentThreshold[1]) + (res$test2 > stringentThreshold[2]) + (res$test3 > stringentThreshold[3]) + (res$test4 < stringentThreshold[4])
  
  return(res)
}


# ##################################
# epigenetic silencing on a per cpg basis as specified in tcga ovarian paper (doi:10.1038/nature10166)
# set by default to exclude NAs since some patients have missing values for methylation
# looking at multiple cpgs per gene can be done by calling this multiple times
#   - geneExpr - a vector of gene level expression for tumor samples only
#   - cpgMethMat - a vector of methylation beta values for tumor samples only
epigeneticSilencingKmeansClusters <- function(geneExpr, cpgMeth){
  
  standMeth <- (cpgMeth - mean(cpgMeth, na.rm=T))/sd(cpgMeth, na.rm=T)
  standExpr <- (geneExpr - mean(geneExpr, na.rm=T))/sd(geneExpr, na.rm=T)
  
  myMat <- cbind(standMeth, standExpr)
  rownames(myMat) <- colnames(methTumor)
  myMat <- na.exclude(myMat)
  
  myK <- kmeans(myMat, centers=2)
  silenced <- ifelse(abs(diff(myK$centers[1, ])) > abs(diff(myK$centers[2, ])), 1, 2)
  silenced <- ifelse(myK$cluster == silenced, 1, 0)
  return(silenced)
}

