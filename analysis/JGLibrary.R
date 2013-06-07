groupMatch <- function(...){
  args <- list(...)
  if(length(args) == 1 && is.list(args[[1]])){
    args <- args[[1]]
  }
  if(length(args) < 2){
    warning("Expected more than one argument")
    return(NULL)
  }
  common <- args[[1]]
  for(i in 2:length(args)){
    common <- intersect(common, args[[i]])
  }
  lapply(args, function(x){ match(common, x)})
}

load.gmt.data <- function(gmt.file.path){
  tmp <- readLines(gmt.file.path)
  gsets <- list()
  for(i in 1:length(tmp)){
    t <- strsplit(tmp[i],'\t')[[1]]
    t2 <- unlist(sapply(t[3:length(t)], function(x){ 
      strsplit(x,"///")
    }),use.names=FALSE)
    gsets[[t[1]]] <- gsub(" ","",t2)
  }
  return (gsets)
}


normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}


# func requires a single parameter which operates on the indices from many
oneToManyOperation <- function(one, many, func=NULL,...){
  idxs <- match(many, one)
  N <- length(unique(one))
  stopifnot(N == length(one))
  map <- vector("list", N)
  for(pos in 1:length(many)){
    if(!is.na(idxs[[pos]])) {
      map[[idxs[pos]]] <- c(map[[idxs[pos]]], pos)
    }
  }
  if(is.null(func)) { return (map)}
  else{ return (sapply(map, func, ...)) }
}

combine_probes_2_gene <- function(expr, genes, method="svd"){
  
  if(is.list(genes)) genes <- unlist(genes)
  
  stopifnot(dim(expr)[1] ==  length(genes))
  ugenes <- sort(setdiff(unique(genes),c("",NULL,NA)))
  
  M <- t(oneToManyOperation(ugenes, genes, function(idxs){
    sub.expr <- as.matrix(expr[idxs,])
    if(dim(sub.expr)[2] > 1){
      tmp <- svd(sub.expr - rowMeans(sub.expr))$v[,1]
      tmp.c <- mean(cor(tmp, t(sub.expr)))
      #cat(gene," ", tmp.c, "\n")
      multiplier <- ifelse(tmp.c < 0, -1, 1)
      sub.expr <- tmp * multiplier
    }
    sub.expr
  }))
  rownames(M) <- ugenes
  colnames(M) <- colnames(expr)
  M
}

fastLMAssociationTest <- function(y, model_matrix){
  p <- ncol(model_matrix)
  n <- nrow(model_matrix)
  
  QR <- qr(model_matrix)
  sig <- sqrt(sum(qr.resid(QR,y)^2) / (n-p+1))
  R <- qr.R(QR)
  xx <- diag(solve(t(R) %*% R))
  
  se <- sig * sqrt(xx)
  2 * pt(abs(qr.coef(QR, y) / se), n-p,lower.tail=FALSE)
}

# performs a fast nested model comparison (3 times faster than standard lm)
fastLMNestedModelComparison <- function(y, model_matrix_1, model_matrix_2){
  
  p1 <- dim(model_matrix_1)[2]
  p2 <- dim(model_matrix_2)[2]
  n1 <- dim(model_matrix_1)[1]
  n2 <- dim(model_matrix_2)[1]
  
  if(n1 != n2){
    stop("Number of samples must be equal.\n")
  }
  
  RSS_1 <- sum(qr.resid(qr(model_matrix_1),y)^2)
  RSS_2 <- sum(qr.resid(qr(model_matrix_2),y)^2)
  
  if(p2 > p1){
    F <- ((RSS_1 - RSS_2) / (p1-p2)) / (RSS_2 / n1-p2)
    return (pf(F, p2-p1,n1-p2,lower.tail=FALSE))
  }else{
    F <- ((RSS_2 - RSS_1) / (p2-p1)) / (RSS_1 / n1-p1)
    return (pf(F, p1-p2,n1-p1,lower.tail=FALSE))
  }
}



