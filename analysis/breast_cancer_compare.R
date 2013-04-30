###############
## compare joe gray, sanger, ccle breast cancer drug response data

heider <- read.table("data~/heider_pnas_breastcancer.txt",sep="\t",header=T,as.is=T,row.names=1)

ccle.breast <- ccle.pharma[getTissueType(rownames(ccle.pharma)) == "BREAST",]
sanger.breast <- sanger.pharma[getTissueType(rownames(sanger.pharma)) == "BREAST",]
cl.ccle <- gsub("(.*?)_.*","\\1",rownames(ccle.breast))
cl.sanger <- gsub("(.*?)_.*","\\1",rownames(sanger.breast))

common_sanger_ccle <- intersect(cl.ccle, cl.sanger)
common_heider_ccle <- intersect(rownames(tbl), cl.ccle)
common_heider_sanger <- intersect(rownames(tbl), cl.sanger)


drug1 <- sanger.pharma[idxs1,"ERLOTINIB"]
drug2 <- tbl[idxs2,"Erlotinib"]

drug.pairs <-  rbind(c("AZD6244","SELUMETINIB","AZD6244"),
                     c("SORAFENIB","SORAFENIB","Sorafenib"),
                     c("PACLITAXEL","PACLITAXEL","Paclitaxel"),
                     c("ERLOTINIB","ERLOTONIB","Erlotinib"),
                     c("NUTLIN-3A","NUTLIN3","Nutlin.3a"),
                     c("LAPATINIB","LAPATINIB","Lapatinib"),
                     c("X17-AAG","DEMETHOXYGELDANAMYCIN","X17.AAG"))

pdf("foo.pdf",width=12,height=20)
par(mfrow=c(nrow(drug.pairs),3))
for(i in 1:nrow(drug.pairs)){
  for(j in 1:3){
    if(j == 1){
      idxs1 <- match(common_sanger_ccle, cl.ccle)
      idxs2 <- match(common_sanger_ccle, cl.sanger)
      drug1 = ccle.pharma[idxs1,drug.pairs[i,2]]
      drug2 = -1 * sanger.pharma[idxs2,drug.pairs[i,1]]
      
    }else if(j==2){
      idxs1 <- match(common_heider_ccle, cl.ccle)
      idxs2 <- match(common_heider_ccle, rownames(heider))
      drug1 = ccle.pharma[idxs1,drug.pairs[i,2]]
      drug2 = heider[idxs2,drug.pairs[i,3]]
    }else{
      idxs1 <- match(common_heider_sanger, cl.sanger)
      idxs2 <- match(common_heider_sanger, rownames(heider))
      drug1 = -1 * sanger.pharma[idxs1,drug.pairs[i,1]]
      drug2 = heider[idxs2,drug.pairs[i,3]]
    }
    N = sum(!is.na(drug1) & !is.na(drug2))
    plot(drug1, drug2,xlab="",ylab="",main=drug.pairs[i,1])
    mtext(paste("R=",format(cor(drug1,drug2,use="pairwise.complete.obs",method="spearman"),digits=2),", N=",N,sep=""))
  }
}
dev.off()