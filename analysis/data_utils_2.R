################

SYN_CCLE_EXPR_ID = "48344"
SYN_CCLE_DRUGRESPONSE_ID = "48359"

synapseLoginFlag="SYN_LOGIN_FLAG"
library("synapseClient")
library("Biobase")
library("affxparser")
library("org.Hs.eg.db")


getSangerMut <- function(){
  tbl <- read.table("data~/sanger_cellline_mutation_data.txt",sep="\t",header=T, as.is=T, quote="",fill=T)
  foo <- t(apply(tbl, 1, function(x){
    apply(cbind(x),1, function(tmp){
      strsplit(tmp,"::")[[1]][1]
    })
  }))
  d <- dim(foo)[2]
  foo <- foo[, c(-2:-5, -(d-c(0:2)))]
  mask <- apply(foo, 2, function(x) all(x=="na"))
  muts <- foo[, !mask]
  muts
}


getCCLEOriginal <- function(){
  loadEntity(getEntity(SYN_CCLE_EXPR_ID))$objects$exprSet
}

getSangerOriginal <- function(){
  loadEntity("syn210931")$objects$eSet_expr
}

getSangerDrugAUC <- function(){
  tbl <- read.table("data~/gdsc_manova_input_w2.txt",sep="\t",header=TRUE,as.is=TRUE,quote="")
  rownames(tbl) <- gsub("-","",toupper(tbl$Cell.Line))
  idxs <- which(grepl("_AUC$", colnames(tbl)))
  auc.tbl <- tbl[,idxs]
  return(tbl)
}

getCharlesSangerDrug <- function(){
  sanger_all <- loadEntity("syn1690438")
  sanger_all <- sanger_all$objects$sanger_data
  assign(x=names(sanger_all)[1],sanger_all[[1]])
  assign(x=names(sanger_all)[2],sanger_all[[2]])
  assign(x=names(sanger_all)[3],sanger_all[[3]])
  assign(x=names(sanger_all)[4],sanger_all[[4]])
  assign(x=names(sanger_all)[5],sanger_all[[5]])
  assign(x=names(sanger_all)[6],sanger_all[[6]])
  a <- rownames(sanger_drug)
  b <- sapply(strsplit(x=colnames(sanger_exp),split="_"),function(x){x[[1]]})
  tmp <- match(a,b)
  sanger_drug$cell_new_id <- colnames(sanger_exp)[tmp]
  sanger_drug <- sanger_drug[!is.na(sanger_drug$cell_new_id),]
  rownames(sanger_drug) <- sanger_drug$cell_new_id
  rm(a,b,tmp)
  sanger_drug$cell_new_id <- NULL
  colnames(sanger_drug) <- toupper(sub(pattern="_IC_50",replacement="",x=colnames(sanger_drug)))
  colnames(sanger_drug) <- sub(pattern=".",replacement="-",x=colnames(sanger_drug),fixed=TRUE)
  tmp <- rownames(sanger_drug)
  sanger_drug <- apply(sanger_drug,2,function(x){as.numeric(sub(pattern=",",replacement=".",x=x,fixed=TRUE))})
  rownames(sanger_drug) <- tmp
  rm(tmp)
  return (sanger_drug)
}

# return pharma and expression as single eset
getCCLE_MetaGenomics <- function(){
  makeExprDrugEset(getCCLEPharma_MetaGenomics(), getCCLEExpr_MetaGenomics())  
}

# return pharma and expression as single eset
getSanger_MetaGenomics <- function(){
  makeExprDrugEset(getSangerPharma_MetaGenomics(), getSangerExpr_MetaGenomics())  
}

getCCLEHybridCapture <- function(){
  ccle_maf <- read.table("data~/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf", header=TRUE,sep="\t",quote="",as.is=TRUE)
  ccle_maf
}

getCCLEHybridCaptureAsMatrix <- function(){
  ccle_maf <- getCCLEHybridCapture()
  genes <- unique(ccle_maf$Hugo_Symbol)
  ids <- unique(ccle_maf$Tumor_Sample_Barcode)
  M <- matrix(0, nrow=length(genes),ncol=length(ids), dimnames=list(genes,ids))
  idxs <- oneToManyOperation(genes, ccle_maf$Hugo_Symbol)
  for(i in 1:length(idxs)){
    M[i,unique(ccle_maf[idxs[[i]],"Tumor_Sample_Barcode"])] <- 1
  }
  return (M)
}

getCCLEPharma_MetaGenomics <- function(){
  getPharma_MetaGenomics('syn1412540')
}

getSangerPharma_MetaGenomics <- function(){
  #getPharma_MetaGenomics('syn1412515')
  sanger.drug <-  pData(loadEntity("syn220680")$objects$adf_drug)
  cellLines <- loadEntity('syn1417611')$objects$cellLines
  m <- matchColumns(data.frame(NAME=rownames(sanger.drug)), cellLines)
  mask <- m != 'notAvailable'
  sanger.drug.m <- sanger.drug[mask,]
  rownames(sanger.drug.m) <- m[mask]
  sanger.drug.m
}

getPharma_MetaGenomics <- function(synapseId){
  cellLines <- loadEntity('syn1417611')$objects$cellLines
  
  # Load the information on the pharmacologic profiling data
  pharma <- loadEntity(synapseId)$objects$pharma
  
  compounds <- read.delim("resources/compounds.txt",sep="\t",header=FALSE, stringsAsFactors=FALSE)
  compounds[which(is.na(compounds[,2])),2] <- compounds[which(is.na(compounds[,2])),1]
  compounds[ compounds[,2] == "",2] <- compounds[ compounds[,2] == "",1]
  compounds[,1] <- toupper(compounds[,1])
  compounds[,2] <- toupper(compounds[,2])
  sym2name <- c(split(compounds[,2], compounds[,1]), split(compounds[,2], compounds[,2]))
  
  inCommon <- intersect(pharma$DRUG, names(sym2name))
  for(i in 1:length(inCommon)){
    pharma$DRUG[which(pharma$DRUG == inCommon[i])] <- sym2name[[inCommon[i]]]
  }
  
  map <- matchColumns(pharma, cellLines)
  pharma$NAME <- map
  
  mat <- matrix(NA, nr=length(unique(pharma$NAME)), nc=length(unique(pharma$DRUG)))
  rownames(mat) <- unique(pharma$NAME); colnames(mat) <- unique(pharma$DRUG)
  
  for(i in 1:nrow(pharma)){
    mat[pharma$NAME[i], pharma$DRUG[i]] <- pharma$activityArea[i]
  }
  
  mat
}

getSangerExpr_MetaGenomics <- function(){
  getExpr(mapId='syn1417190',dataId='syn427896') 
}


getCCLEExpr_MetaGenomics <- function(){
  getExpr(mapId='syn1417166',dataId='syn425002')
}

getExpr <- function(mapId, dataId){
  cellLines <- loadEntity('syn1417611')$objects$cellLines
  
  # Load the various data matrices
  mapDat <- loadEntity(mapId)$objects$mapDat
  data <- loadEntity(dataId)
  dat <- exprs(data$objects$eset)
  
  # Find the columns in the mapDat objects corresponding to the column names
  mapDat <- mapDat[match(colnames(dat), mapDat$technologyID),]
  
  map <- matchColumns(mapDat, cellLines)
  
  expr <- averageReplicates(dat, map)
  
  entrez.ids <- gsub("(\\d*)_mt","\\1", rownames(expr))
  genes <- unlist(mget(entrez.ids, org.Hs.egSYMBOL, ifnotfound=NA))
  
  tmp <- combine_probes_2_gene(expr, genes)
  eset <- new("ExpressionSet",expr=tmp)
  eset
}


averageReplicates <- function(dat,map){
  ids <- which(map != "notAvailable")
  dat <- dat[,ids]
  map <- map[ids]
  map2row <- split(1:ncol(dat), map)
  tmp <- matrix(NA, nr=nrow(dat),nc=length(map2row))
  colnames(tmp) <- names(map2row)
  rownames(tmp) <- rownames(dat)
  for(i in 1:length(map2row)){
    cat("\r",i)
    if(length(map2row[[i]]) == 1){
      tmp[,i] <- dat[,map2row[[i]]]
    }else{
      tmp[,i] <- rowMeans(dat[,map2row[[i]]])
    }
  }
  tmp
}


matchColumns <- function(dat, cellLines=NULL){
  if(is.null(cellLines)){ cellLines <- loadEntity('syn1417611')$objects$cellLines}
  ids1 <- match(toupper(dat$cellLineName), toupper(cellLines$Cell.line.name))
  ids2 <- match(toupper(dat$cellLineName), toupper(cellLines$Simplified.cell.name))
  ids3 <- match(toupper(dat$NAME), toupper(cellLines$Cell.line.name))
  ids4 <- match(toupper(dat$NAME), toupper(cellLines$Simplified.cell.name))
  mat <- cbind(ids1,ids2,ids3,ids4)
  newNames <- apply(mat, 1, function(x){ 
    if(sum(!is.na(x)) == 0){ 
      return("notAvailable")
    }else{ 
      cellLines$TIP_name[x[which(!is.na(x))[1]]]
    }
  })
  
  newNames
}
combine_probes_2_gene <- function(expr, genes, method="svd"){
  
  if(is.list(genes)) genes <- unlist(genes)
  
  stopifnot(dim(expr)[1] ==  length(genes))
  ugenes <- unique(genes)
  ugenes <- sort(ugenes[!is.na(ugenes)])
  M <- matrix(NaN, ncol=dim(expr)[2],nrow=length(ugenes),
              dimnames=list(ugenes, colnames(expr)))
  
  for(gene in ugenes){
    sub.expr <- as.matrix(expr[which(genes == gene),])
    if(dim(sub.expr)[2] == 1){
      M[gene,] <- sub.expr
    }else{
      tmp <- svd(sub.expr - rowMeans(sub.expr))$v[,1]
      tmp.c <- mean(cor(tmp, t(sub.expr)))
      #cat(gene," ", tmp.c, "\n")
      multiplier <- ifelse(tmp.c < 0, -1, 1)
      M[gene,] <- tmp * multiplier
    }
  }
  M
}

getTissueType <- function(dat){
  cn <- dat
  if(is.matrix(dat)){
    cn <- colnames(dat)
  }
  if(!is.vector(cn)) stop("Expected vector")
  gsub(".*?_(.*)","\\1", cn, perl=TRUE)
}

makeExprDrugEset <- function(pharmaDat, exprDat){
  idxs <- match(rownames(pharmaDat), sampleNames(exprDat))
  pharmaDat.m <- as.data.frame(pharmaDat[!is.na(idxs),])
  exprDat.m <- exprDat[, na.omit(idxs)]
  eset <- new("ExpressionSet", expr=exprs(exprDat.m), phenoData=new("AnnotatedDataFrame",data=pharmaDat.m))
  eset
}

