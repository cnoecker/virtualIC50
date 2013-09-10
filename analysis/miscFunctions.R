library(e1071) # for hamming distance

compute_residue_locations <- function(aachanges){
  lapply(aachanges, function(x){
    if(x == "") { return (0)}
    sapply(strsplit(x,";")[[1]], function(y){
      as.numeric(gsub("p\\.\\w(\\d+).*","\\1",y))
    })
  })
}

egfr_domains <- function(aachanges){
  domains <- sapply(aachanges, function(x){
    if(x == "") { return ("WT")}
    tmp <- sapply(strsplit(x,";")[[1]], function(y){
      pos <- as.numeric(gsub("p\\.\\w(\\d+).*","\\1",y))
      cat(pos,"\n")
      if(pos < 646){ return ("extraceullar")}
      else if(pos >= 646 & pos < 669) { return ("helical")}
      else if(pos >= 669) { return ("cytoplasmic")}
    })
    return (paste(sort(tmp),collapse=";"))
  })
  return (domains)
}

pik3ca_domains <- function(aachanges){
  domains <- sapply(aachanges, function(x){
    if(x == "") { return ("WT")}
    tmp <- sapply(strsplit(x,";")[[1]], function(y){
      pos <- as.numeric(gsub("p\\.\\w(\\d+).*","\\1",y))
      #cat(pos,"\n")
      if(pos < 16){ return ("ABD-")}
      else if(pos >= 16 & pos < 106) { return ("ABD")}
      else if(pos >= 106 & pos < 187) { return ("ABD-RBD")}
      else if(pos >= 187 & pos < 290) { return ("RBD")}
      else if(pos >= 290 & pos < 330) { return ("RBD-C2")}
      else if(pos >= 330 & pos < 487) { return ("C2")}
      else if(pos >= 487 & pos < 517) { return ("C2-helical")}
      else if(pos >= 517 & pos < 695) { return ("helical")}
      else if(pos >= 695 & pos < 797) { return ("helical-kinase")}
      else if(pos >= 797) { return ("kinase")}
    })
    return (paste(sort(tmp),collapse=";"))
  })
  return (domains)
}

## expects a logical or numeric (0,1) features by samples matrix
combine.features.by.hammingdistance <- function(M,mincountPerRow=3,hammingDistThreshold=1){ 

  Mm <- M[rowSums(M) >= mincountPerRow, ]
  hd <- hamming.distance(Mm)
  hd[lower.tri(hd,diag=TRUE)] <- 100 # mask out lower triangle
  idxs <- which(hd <= hammingDistThreshold, arr.ind=TRUE)
  if(length(idxs) == 0){
    return (Mm)
  }
  groupings <- list()
  for(i in 1:nrow(idxs)){
    idxpair <- idxs[i,]
    found = FALSE
    if(length(groupings) > 0){
      for(j in 1:length(groupings)){
        if(any(idxpair %in% groupings[[j]])){
          groupings[[j]] <- union(groupings[[j]], idxpair)
          found = TRUE
        }
      }
    }
    if(!found){
      groupings[[length(groupings)+1]] <- idxpair
    }
  }
  
  rmIdxs <- unique(do.call("c", groupings))
  metaM <- do.call("rbind",lapply(groupings, function(x) { as.numeric(colSums(Mm[x,]) > 0) }))
  metaGenes <- sapply(groupings, function(x){ paste(sort(rownames(hd)[x]),collapse="::") })
  rownames(metaM) <- metaGenes
  tmp <- rbind(Mm[-rmIdxs,], metaM)
  return (tmp)
}
