##############
# semi-supervised implementation of elastic net: Mark Culp "On The Semi-Supervised Joint Trained Elastic Net"

bootstrappedJointTrainingElasticNet <- function(Xl, Xu, y, Xv=list(), alphas=.1, gamma=1,lambda.min=TRUE,scale=TRUE,numbootstraps=10,ncores=0){
  require(glmnet)
  
  Nu <- nrow(Xu)
  
  if(scale){
    Xl <- scale(Xl)
    Xu <- scale(Xu)
    Xv <- lapply(Xv, function(x) scale(x))
  }
  y_c <- c(y - mean(y), rep(0, Nu))
  
  
  if(gamma == 1){
    Xu_g <- Xu
  }else{
    pc <- svd(Xu)
    tmp <- t(pc$u) %*% Xu
    P <- 1 / sqrt(pc$d^2 / gamma + 1)
    Xu_g <- diag(P) %*% tmp 
  }
  N <- nrow(Xl)
  fits <- mclapply(1:numbootstraps, function(i){
    idxs <- sample(N,replace=TRUE)
    X <- rbind(Xl[idxs,], Xu_g)
    y_c <- c(y[idxs] - mean(y[idxs]), rep(0, Nu))
    cv <- cv.glmnet(X,y_c,alpha=alpha,nfolds=5,standardize=FALSE)
    lambda <- ifelse(lambda.min, cv$lambda.min, cv$lambda.1se)
    fit <- glmnet(X,y_c,alpha=alpha,lambda=lambda,standardize=FALSE)
    fit 
  },mc.cores=ncores,mc.set.seed=TRUE,mc.preschedule=FALSE)
  
  yhats <- lapply(Xv, function(x){ 
    preds <- sapply(fits, function(fit) predict(fit, x) )
    rowMeans(preds)
  })
  
  return(yhats)
}

jointTrainingElasticNet <- function(Xl, Xu, y, Xv=list(), alphas=.1, gammas=0,lambda.min=TRUE,scale=TRUE,compare.regularEN=FALSE){
  require(glmnet)

  Nu <- nrow(Xu)
  
  if(scale){
    Xl <- scale(Xl)
    Xu <- scale(Xu)
    Xv <- lapply(Xv, function(x) scale(x))
  }
  y_c <- c(y - mean(y), rep(0, Nu))
  #browser()
  pc <- svd(Xu)
  tmp <- t(pc$u) %*% Xu
  
  R <- vector(mode = "list", length = length(gammas) * length(alphas) + 1)
  
  for(i in seq_along(gammas)){
    gamma <- gammas[i]
    if(gamma == 0){
      P <- rep(0, Nu)
    }else{
      P <- 1 / sqrt(pc$d^2 / gamma + 1)
    }
    Xu_g <- diag(P) %*% tmp 
    X <- rbind(Xl, Xu_g)
    for(j in seq_along(alphas)){
      alpha <- alphas[j]
      cv <- cv.glmnet(X,y_c,alpha=alpha,nfolds=5,standardize=FALSE)
      lambda <- ifelse(lambda.min, cv$lambda.min, cv$lambda.1se)
      fit <- glmnet(X,y_c,alpha=alpha,lambda=lambda,standardize=FALSE)
      yhats <- lapply(Xv, function(x){ predict(fit, x) })
      idx <- length(gammas) * (j-1) + i
      R[[idx]] <- list(alpha=alpha,gamma=gamma,lambda=lambda,yhats=yhats)
    }
  }
  
  if(compare.regularEN){
    cv <- cv.glmnet(Xl,y-mean(y),alpha=.1,nfolds=5)
    lambda <- ifelse(lambda.min, cv$lambda.min, cv$lambda.1se)
    fit <- glmnet(Xl,y-mean(y),alpha=.1,lambda=lambda)
    R[[length(R)]] <- list(yhats=lapply(Xv, function(x){ predict(fit, x) }))
  }
  
  R
}


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
