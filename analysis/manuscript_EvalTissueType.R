source("analysis/JGLibrary.R")
source("analysis/data_utils_2.R")
source("analysis/celline_2_tcga_pipeline.R")
source("analysis/miscFunctions.R")
source("analysis/cellline_2_tcga_mutationImporter.R")
library(synapseClient)
source("code/ssVDSFunctions.R")

ncores=10

if(!exists("ccle",.GlobalEnv)){
  #synapseLogin("justin.guinney@sagebase.org",'marley')
  #ccle <- getCCLE_MetaGenomics()
  load("~/data/ccle.rda")
  load("~/data/sanger.rda")
}

if(!exists("tcgaList",.GlobalEnv)){
  env <- new.env()
  load("~/data/TCGA_ds_ver4.rda", envir=env)
  tcgaList <- as.list(env)
}

ccle_drugs <- c("ERLOTONIB","SELUMETINIB","LAPATINIB","PACLITAXEL","SORAFENIB","TKI258","PLX4720")
sanger_drugs <- c("Cisplatin","Temsirolimus","AZD-2281","Sunitinib","Doxorubicin")
all_drugs <- c(ccle_drugs, sanger_drugs)
is.sanger <- c(rep(FALSE, length(ccle_drugs)),rep(TRUE, length(sanger_drugs)))
diseases <- c("luad","lusc","skcm","ov","crc","brca","brca.3neg","ucec","kirc","stad")
cldiseases <- c("LUAD","LUSC","SKIN","OVARY","LARGE_INTESTINE","BREAST","BREAST_3neg","ENDOMETRIUM", "KIDNEY","STOMACH")


## prepare ccle data
ccleNoLymph <- ccle[,!(getTissueType(sampleNames(ccle)) %in% c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))]
load("/home/cferte/resources/cell_line_infos.rda")
ccleinfo <- cells.info$ccle.info
tt <- getTissueType(sampleNames(ccleNoLymph))
histology <- ccleinfo[sampleNames(ccleNoLymph),"Hist.Subtype1"]
clnames <- gsub("(.*?)_.*","\\1", sampleNames(ccleNoLymph), perl=TRUE)
luad.mask <- tt == "LUNG" & histology == "adenocarcinoma"
lusc.mask <- tt == "LUNG" & histology == "squamous_cell_carcinoma"
lbls <- sampleNames(ccleNoLymph)
lbls[luad.mask] <- paste(clnames[luad.mask],"_LUAD",sep="")
lbls[lusc.mask] <- paste(clnames[lusc.mask],"_LUSC",sep="")
sampleNames(ccleNoLymph) <- lbls

## prepare sanger data
sangerNoLymph <- sanger[,!(getTissueType(sampleNames(sanger)) %in% c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))]
sangerinfo <- cells.info$sanger.info

clnames <- gsub("(.*?)_.*","\\1", sampleNames(sangerNoLymph), perl=TRUE)
idxs <- match(clnames,toupper(gsub("-","", rownames(sangerinfo))))
hist <- sangerinfo[idxs,"Tissue"]

luad.mask <- grepl("lung: NSCLC: adenocarcinoma",hist)
lusc.mask <- grepl("lung: NSCLC: squamous_cell_carcinoma",hist)
lbls <- sampleNames(sangerNoLymph)
lbls[luad.mask] <- paste(clnames[luad.mask],"_LUAD",sep="")
lbls[lusc.mask] <- paste(clnames[lusc.mask],"_LUSC",sep="")
sampleNames(sangerNoLymph) <- lbls


tcgaList <- tcgaList[diseases]
for(tt in names(tcgaList)){ tcgaList[[tt]]$fmat <- build_feature_matrix(tcgaList[[tt]],with.rppa=FALSE) }
tcgaGeneExpr <- lapply(tcgaList, function(x) x$geneExpr)

outdir="./results/paper_2013_08_20/"

for(drug_idx in seq_along(all_drugs)){
  drug <- all_drugs[drug_idx]
  sanger_flag <- is.sanger[drug_idx]
  eset <- if(sanger_flag){ sangerNoLymph } else { ccleNoLymph }
  vds <- virtual_ic50(eset, drug, tcgaGeneExpr, seed=2013, reverseDrug=sanger_flag,perf.eval=TRUE,num.bootstraps=20,ncores=ncores,alpha=.1)
  save(vds, file=paste(outdir,"rda/",drug,".rda",sep=""))
  
  for(disease in names(vds$yhat)){
    yhat <- vds$yhat[[disease]]$y_mean
    #disease <- names(tcgaList)[tcga_idx]
    fmat <- find_drug_features(yhat,tcgaList[[disease]]$fmat,rep(TRUE, nrow(tcgaList[[disease]]$fmat)), beta_threshold=10^-3,ncores=ncores,num.bootstraps=100, outlier.sd=8)
    pdf(paste(outdir,"bubblePlots/",drug,"_",disease,".pdf",sep=""),width=10,height=6,useDingbats=F)
    plot_features_3(fmat,paste("TCGA ",disease,sep=""),top=35,text.cex=.7)
    dev.off()
    save(fmat, file=paste(outdir,"rda/",drug,"_",disease,".rda",sep=""))
  }
}

################################################
######## Generate performance plots

#all_rda <- dir(paste(outdir,"rda",sep=""))
#disease_rda <- dir(paste(outdir,"rda",sep=""),pattern="_")
#cell_rda <- setdiff(all_rda, disease_rda)


for(drug in all_drugs){
  drugpath <- paste(outdir,"rda/",drug,".rda",sep="")
  if(!file.exists(drugpath)){ next; }
  
  
  env <- new.env()
  load(drugpath,envir=env)
  cltypes <- names(env$vds$perf)
  
  all.tt <- sort(union(cldiseases, cltypes))
  clAucs <- list()
  clAucs[all.tt] <- 0
  tcgaAucs <- list()
  tcgaAucs[all.tt] <- 0
  
  aucs <- sapply(env$vds$perf, function(x) x$auc)
  clAucs[names(aucs)] <- aucs
  #for(cl in all.tt){ clAucs[[cl]] <- env$vds$perf[[cl]]$auc }
  #rho <- sapply(cltypes, function(x)env$vds$perf[[x]]$rho$estimate)
  
  for(idx in seq_along(diseases)){
    tcgaName <- diseases[idx]
    disease <- cldiseases[idx]
    drugdiseasepath <- paste(outdir,"rda/",drug,"_",tcgaName,".rda",sep="")
    if(!file.exists(drugdiseasepath)){ next; }
    env <- new.env()
    load(drugdiseasepath,envir=env)
    tcgaAucs[[disease]] <- env$fmat$metric["auc"]
  }
  
  b2r <- heat.colors(20)[cut(1-unlist(clAucs), seq(0,1,by=.05),include.lowest=T)]
  pdf(paste(outdir,"perfPlots/",drug,".pdf",sep=""),width=10,height=8)
  par(mfrow=c(1,2))
  par(oma=c(3,10,1,1))
  barplot(unlist(clAucs), col=b2r, las=2,horiz=TRUE,xlim=c(0,1),xlab="AUC");
  abline(v=c(.2,.4, .6 ,.8),lty=2)
  abline(v=.5,lwd=2)
  
  b2r <- heat.colors(20)[cut(1-unlist(tcgaAucs), seq(0,1,by=.05),include.lowest=T)]
  #par(oma=c(1,10,1,1))
  barplot(unlist(tcgaAucs), col=b2r, las=2,horiz=TRUE,xlim=c(0,1),xlab="AUC",yaxt="n");
  abline(v=c(.2,.4,.6 ,.8),lty=2)
  abline(v=.5,lwd=2)
  dev.off()
}

############################
#### generate pathways based on enriched kegg pathways
library("KEGG.db")
library(org.Hs.eg.db)
library(pathview)

keggPwys <- as.list(KEGGPATHID2EXTID)
keggPwys <- keggPwys[grepl("hsa", names(keggPwys))]
lpwys <- sapply(keggPwys, length)
keggPwys <- keggPwys[lpwys >= 10 & lpwys < 500]
keggPwyNames <- unlist(mget(gsub("hsa","",names(keggPwys)), KEGGPATHID2NAME))
idxs  <- grep("signaling|excision|cancer|replication|autophagy|interaction|Spliceosome|cycle|Apoptosis", keggPwyNames)

cancerPwys <- keggPwys[idxs]
cancerPwyNames <- keggPwyNames[idxs]


rda.files <- list.files(paste(outdir,"rda/",sep=""),pattern="*.rda")
rda.files <- rda.files[grepl("_", rda.files)]

oncogenes <- getOncoGenes()
oncogenesAsEntrez <- na.omit(mget( oncogenes, org.Hs.egSYMBOL2EG,ifnotfound=NA))

pwydir <- paste(outdir,"kegg",sep="")
dir.create(pwydir)

top.traits <- 25
pval.adj.threshold = .05

for(f in rda.files){
  drugdiseasepath <- paste(outdir,"rda/",f,sep="")
  env <- new.env()
  load(drugdiseasepath,envir=env)
  traits <- rownames(env$fmat$df)[1:top.traits]

  IS <- (env$fmat$df$freqCounts * sign(env$fmat$df$effect))[1:top.traits]
  
  genesL <- lapply(traits, function(x){ 
    sapply(strsplit(x,"::")[[1]], function(y) strsplit(y,"_")[[1]][1]) } )
  IS.v <- unlist(sapply(1:top.traits, function(i) {
    rep(IS[i], length(genesL[[i]]))
  } ))
  names(IS.v) <- unlist(genesL)
  entrez.ids <- unlist(mget(unlist(genesL), org.Hs.egSYMBOL2EG,ifnotfound=NA))
  pvals <- sapply(cancerPwys, function(x){ 
    fisher.test(factor(oncogenesAsEntrez %in% entrez.ids,levels=c("FALSE","TRUE")),
                factor(oncogenesAsEntrez %in% x,levels=c("FALSE","TRUE")),alternative="greater")$p.value})
  
  sigPwys <- names(which(p.adjust(pvals,method="BH") < pval.adj.threshold))
  for(pwy in sigPwys){
    current.dir <- getwd()
    setwd(pwydir)
    pv.out=pathview(gene.data=IS.v,pathway.id=pwy,
                    out.suffix=f,kegg.native=T,
                    gene.idtype="SYMBOL",gene.annotpkg="org.Hs.eg.db",
                    both.dirs=list(gene=T,cpd=T),same.layer=FALSE)
    #pathview(gene.data=v,pathway.id="05200",gene.idtype="SYMBOL",gene.annotpkg="org.Hs.eg.db") 
    setwd(current.dir)
  }
}

###############
### logic regression pipeline
library(LogicReg)

run.logicreg <- function(rdafile,nleaves=c(1,10),min.freq=.05,top.features=25){
  env <- new.env()
  load(rdafile,envir=env)
  top.features <- rownames(env$fmat$df)[1:top.features]
  vds <- env$fmat$vds
  A <- env$fmat$dataMatrix
  idxs <- match(names(vds), colnames(A))
  vds.m <- vds[!is.na(idxs)]
  A.m <- A[, na.omit(idxs)]
  freq <- rowSums(A.m) / ncol(A.m)
  dm <- A.m[which(rownames(A.m) %in% top.features & freq > min.freq),,drop=FALSE]
  if(nrow(dm) > 1){
    myanneal <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 0)
    fit <- logreg(vds.m, t(dm), select=2,nleaves=nleaves, ntree=1,anneal.control = myanneal)
    cvfit <- logreg(vds.m, t(dm), select=3,nleaves=nleaves, ntree=1,anneal.control = myanneal)
    
    scores <- cvfit$cvscores$test.ave[cvfit$kfold * seq(cvfit$nleaves[1],cvfit$nleaves[2])]
    best.model.idx <- sort(which(scores == min(scores)))[1]
    out <- capture.output(print(fit, nms=rownames(dm)))
    
    return (list(fit=fit,cvfit=cvfit,all.models=out,best.model=out[best.model.idx*2]))
  }else{
    return (NULL)
  }
}

rda.files <- list.files(paste(outdir,"rda/",sep=""),pattern="*.rda")
rda.files <- rda.files[grepl("_", rda.files)]
for(f in rda.files){
  r <- run.logicreg(paste(outdir,"rda/",f,sep=""))
  if(!is.null(r)){
    pdf(paste(outdir,"logreg/",f,".pdf",sep=""),width=9,height=5)
    op <- par(
      oma=c(0,0,.5,0),# Room for the title and legend
      mfrow=c(1,2)
    )
    plot(r$fit)
    plot(r$cvfit)
    par(op)
    mtext(r$best.model, line=2, font=2, cex=.7)
    # plot all models
    idxs <- seq(2, 20,by=2)
    par(mar=c(0,2,0,0)); 
    plot(1:10, type="n", axes=F, xlab="", ylab=""); text(rep(1, 10), rev(1:10), labels=r$all.models[idxs],cex=.5,adj=0)
    axis(side=2, at=c(1:10), labels=as.character(10:1),las=1)
    dev.off()
  }
}

###############
### meta-pathway analysis

################
### misc ...

sangerMut <- getSangerMut()
tms <- pData(sangerNoLymph)[,"Temsirolimus"]
names(tms) <- colnames(sangerNoLymph)
tms <- tms[!is.na(tms) & grepl("BREAST",names(tms))]
clNames <- gsub("(.*?)_.*","\\1", names(tms))
idxs <- match(clNames, toupper(gsub("-","", sangerMut[,1])))
breastSangerMut <- sangerMut[idxs,]

sangerCNV <- getSangerCNV()
idxs <- match(clNames, toupper(gsub("-","", colnames(sangerCNV))))
breastSangerCNV <- exprs(sangerCNV[, idxs])
myc_amp <- breastSangerCNV["MYC",]
tp53_del <- breastSangerCNV["TP53",]
akt_del <- breastSangerCNV["AKT1",]

tp53_del_cutoff <- -0.72
myc_amp_cutoff <- 3
akt_del_cutoff <- -0.2

sangerDM <- rbind(
  tp53_del < tp53_del_cutoff,
  breastSangerMut[,"TP53"] != "wt",
  akt_del < akt_del_cutoff, 
  breastSangerMut[,"PIK3CA"] != "wt",
  myc_amp > myc_amp_cutoff)
rownames(sangerDM) <- c("TP53_del","TP53_mut","AKT1_del","PIK3CA_mut","MYC_amp")

p <- predict(fit, msz=7, newbin=t(sangerDM))

f1 <- -tms ~ factor(p[,1], labels=c("resistant","sensitive"))
f2 <- -tms ~ factor(breastSangerMut[,"PIK3CA"] != "wt", labels=c("PIK3CA wt","PIK3CA mut"))
f3 <- -tms ~ factor(breastSangerMut[,"TP53"] != "wt", labels=c("tp53 wt","tp53 mut"))
f4 <- -tms ~ factor(myc_amp > myc_amp_cutoff, labels=c("myc wt","myc amp"))
f5 <- -tms ~ factor(akt_del < akt_del_cutoff, labels=c("akt_wt","akt_del"))
f6 <- -tms ~ factor(tp53_del < tp53_del_cutoff, labels=c("tp53_wt","tp53_del"))

par(mfrow=c(2,3))
boxplot(f1,ylab="-IC50"); mtext(paste("p=",format(wilcox.test(f1)$p.value,digits=2),sep=""))
boxplot(f2,ylab="-IC50"); mtext(paste("p=",format(wilcox.test(f2)$p.value,digits=2),sep=""))
boxplot(f3,ylab="-IC50"); mtext(paste("p=",format(wilcox.test(f3)$p.value,digits=2),sep=""))
boxplot(f4,ylab="-IC50"); mtext(paste("p=",format(wilcox.test(f4)$p.value,digits=2),sep=""))
boxplot(f5,ylab="-IC50"); mtext(paste("p=",format(wilcox.test(f5)$p.value,digits=2),sep=""))
boxplot(f6,ylab="-IC50"); mtext(paste("p=",format(wilcox.test(f6)$p.value,digits=2),sep=""))