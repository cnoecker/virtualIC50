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
