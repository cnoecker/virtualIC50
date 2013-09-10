
plot_data_ellipse <- function(x, y,...){
  x0 <- mean(x)
  y0 <- mean(y)
  a <- sd(x)
  b <- sd(y)
  theta <- seq(0, 2 * pi, .01)
  xl <- x0 + a * cos(theta)
  yl <- y0 + b * sin(theta)
  lines(xl, yl, ...)
}

make_tissue_centroid_plot <- function(eset, tissues, dataset){
  tf <- factor(tissues)
  tl <- levels(tf)
  
  pc <- svd(exprs(eset) - rowMeans(exprs(eset)))
  tmp <- sapply(tl, function(tissue) {
     mask <- tf == tissue
     c(mean(pc$v[mask,1]), mean(pc$v[mask,2]))
  })
  
  col = rainbow(length(tl))
  ascii_chars <- c(65:90, 97:122) 
  pdf(paste("plots/",dataset,"_centroid_pcplot.pdf",sep=""),width=8,height=6)
  par(oma=c(0,0,0,10))
  #plot(t(tmp), xlab="PC1", ylab="PC2",
  #     col=col,cex=1.5,pch=ascii_chars[1:length(tl)])
  lbls <- substring(tl, 1, 3)
  plot(t(tmp)*1.2, xlab="PC1", ylab="PC2",type="n")
  text(t(tmp), labels=lbls,col=col)
  for(i in 1:length(tl)){
    mask <- tl[i] == tissues
    if(sum(mask) > 2){
      plot_data_ellipse(pc$v[mask,1], pc$v[mask,2],col=col[i],lwd=1, lty=2)
    }
  }
  
  par(xpd=NA)
  legend(par("usr")[2], par("usr")[4],legend=tl,
         #pch=ascii_chars[1:length(tl)],
         text.col=rainbow(length(tl)),
         cex=.7)
  par(xpd=FALSE)
  dev.off()
}

make_fan_plot <- function(eset, tissues, dataset, var.iqr.filter=.5){
  v <- apply(exprs(eset), 1, var)
  X <- exprs(eset[v > quantile(v, .5),])
  dist.x <- as.dist(1 - cor(X,method="spearman"))
  hc <- hclust(dist.x,method="ward")
  
  
  #pdf("plots/ccle_cordist_dendrogram.pdf",width=30,height=4)
  #ColorDendrogram(hc, y=as.numeric(factor(ccle_tissues)),
  #                labels=ascii_chars[as.numeric(factor(ccle_tissues))],
  #                branchlength=5)
  pdf(paste("plots/", dataset, "_cordist_fan.pdf",sep=""), width=15,height=15)
  hc$labels <- tissues
  plot(as.phylo(hc),type="fan",cex=.6,
       tip.color=rainbow(length(levels(factor(tissues))))[as.numeric(factor(tissues))])
  dev.off()
}


makeTissuePCPlot <- function(eset, tissues, dataset){
  
  ascii_chars <- c(65:90, 97:122) 
  pc <- svd(exprs(eset) - rowMeans(exprs(eset)))
  
  pdf(paste("plots/",dataset,"_pcplot.pdf",sep=""),width=11,height=10)
  par(mfrow=c(2,2),oma=c(10,1,1,1))
  plot(pc$d^2 / sum(pc$d^2),ylab="% var",pch=19,cex=.5,xlab="eigen index")
  par(mar=c(5,10,4,2))
  barplot(sort(table(factor(tissues))),horiz=T, cex.names=0.5,las=2)
  
  tf <- factor(tissues)
  tl <- levels(tf)
  col = rainbow(length(tl))
  
  plot(pc$v[,1], pc$v[,2], xlab="PC1", ylab="PC2",
       col=col[as.numeric(tf)],cex=.5,pch=ascii_chars[as.numeric(tf)])
  par(xpd=NA)
  legend(par("usr")[2], par("usr")[4],legend=tl,
         pch=ascii_chars[1:length(tl)],
         col=col,
         cex=.7)
  par(xpd=FALSE)
  dev.off()
}
