library(synapseClient)
library(Biobase)

getCCLEMaf <- function(){
  ccle_maf <- read.table("data~/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf", header=TRUE,sep="\t",quote="",as.is=TRUE)
  ccle_maf
}

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

getSangerExpr <- function(){
  getExpr(mapId='syn1417190',dataId='syn427896') 
}

getCClEExpr <- function(){
  getExpr(mapId='syn1417166',dataId='syn425002') 
}

getGSKExpr <- function(){
  getExpr(mapId='syn1417186',dataId='syn1417646')
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
  
  averageExpr <- averageReplicates(dat, map)
  averageExpr
}

#Function to annotate column names with uniform identifier
matchColumns <- function(dat, cellLines){
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
  rownames(mat) <- dat$technologyID
  #   cbind(mat,dat, newNames)
  newNames
}

getTissueType <- function(dat){
  cn <- dat
  if(is.matrix(dat)){
    cn <- colnames(dat)
  }
  if(!is.vector(cn)) stop("Expected vector")
  gsub(".*?_(.*)","\\1", cn, perl=TRUE)
}

getCellLineName <- function(dat){
  cn <- dat
  if(is.matrix(dat)){
    cn <- colnames(dat)
  }
  if(!is.vector(cn)) stop("Expected vector")
  gsub("(.*?)_.*","\\1", cn, perl=TRUE)
}

# Function to average replicates
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

getCCLEOriginal <- function(){
  SYN_CCLE_EXPR_ID = "48344"
  ccle_entity_expr <- loadEntity(getEntity(SYN_CCLE_EXPR_ID))$objects$exprSet
  ccle_entity_expr
}

getSangerOriginal <- function(){
  eset <- loadEntity("syn210931")$objects$eSet_expr
  eset
}


