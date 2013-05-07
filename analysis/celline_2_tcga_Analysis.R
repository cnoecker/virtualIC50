source("analysis/data_utils_2.R")
source("analysis/celline_2_tcga_pipeline.R")
library(synapseClient)

synapseLogin("justin.guinney@sagebase.org",'marley')

sanger <- getSanger_MetaGenomics()
ccle <- getCCLE_MetaGenomics()


###############
## build data sets
brca <- build.tcga.ds(geneExprId="syn1446183",rppaId="syn1571265",gisticId="syn1687590",cbioPrefix="brca",isRNASeq=TRUE)
luad <- build.tcga.ds(geneExprId="syn418003",rppaId="syn464306",gisticId="syn1687610",cbioPrefix="luad",isRNASeq=TRUE)
coad <- build.tcga.ds(geneExprId="syn1446195",rppaId="syn1446043",gisticId="syn1687596",cbioPrefix="coadread",isRNASeq=TRUE)
read <- build.tcga.ds(geneExprId="syn1446274",rppaId="syn1446053",gisticId="syn1687628",cbioPrefix="coadread",isRNASeq=TRUE)
skcm <- build.tcga.ds(geneExprId="~/projects/h3/data~/firehose/gdac.broadinstitute.org_SKCM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2013032600.0.0/SKCM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
                      gisticId="syn1687618",
                      rppaId="~/projects/h3/data~/firehose/gdac.broadinstitute.org_SKCM.RPPA_AnnotateWithGene.Level_3.2013032600.0.0/SKCM.rppa.txt",
                      cbioPrefix="skcm",isRNASeq=TRUE)
blca <- build.tcga.ds(geneExprId="syn417761",rppaId="syn1681031",gisticId="syn1687592",cbioPrefix="syn1571577",isRNASeq=TRUE)
lusc <- build.tcga.ds(geneExprId="syn1446244",rppaId="syn1446049",gisticId="syn1687612",cbioPrefix="lusc",isRNASeq=TRUE)


# merge rectal and colon
idxs1 <- groupMatch(rownames(coad$geneExpr), rownames(read$geneExpr))
idxs2 <- groupMatch(rownames(coad$rppa), rownames(read$rppa))
idxs3 <- groupMatch(rownames(coad$gistic), rownames(read$gistic))
crc <- list(geneExpr=cbind(coad$geneExpr[idxs1[[1]],], read$geneExpr[idxs1[[2]],]),
            rppa=cbind(coad$rppa[idxs2[[1]],], read$rppa[idxs2[[2]],]),
            gistic=cbind(coad$gistic[idxs3[[1]],], read$gistic[idxs3[[2]],]),
            mut=coad$mut)

# save datasets for fast retreival
save(brca, crc, blca, luad, skcm, lusc, file="~/data/TCGA_ds.rda")


# triple neg brca
brca.3neg <- brca
brca.3neg$geneExpr <- brca.3neg$geneExpr[,as.matrix(brca$geneExpr)["ESR1",] < 8 & as.matrix(brca$geneExpr)["PGR",] < 8]
brca.3neg$gistic <- brca.3neg$gistic[, as.matrix(brca$gistic)["ERBB2",] != 2]

###############
## fgfr in bladder
carcinoma_mask <- !(getTissueType(sampleNames(sanger)) %in% c("CENTRAL_NERVOUS_SYSTEM","SKIN","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
y_hat_fgfr <- virtual_ic50(sanger[,carcinoma_mask],"PD-173074",blca,reverseDrug=TRUE)
blca_fgfr_vogel <- find_drug_features (y_hat_fgfr,blca, with.rppa=FALSE,beta_threshold=10^-3,gene.dict="vogelstein")
blca_fgfr_cosmic <- find_drug_features (y_hat_fgfr,blca, with.rppa=FALSE,beta_threshold=10^-3,gene.dict="cosmic")
blca_fgfr_cbio <- find_drug_features (y_hat_fgfr,brca, with.rppa=FALSE,beta_threshold=10^-3,gene.dict="cbio")

pdf("plots/brca_fgfr_vogel.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_fgfr_vogel, "PD-173074 (FGFR)","BRCA",text.cex=.9)
dev.off()
pdf("plots/brca_fgfr_cosmic.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_fgfr_cosmic, "PD-173074 (FGFR)","BRCA",text.cex=.9)
dev.off()
pdf("plots/brca_fgfr_cbio.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_fgfr_cbio, "PD-173074 (FGFR)","BRCA",text.cex=.9)
dev.off()


##########################
# BRAF inhib in melanoma
#melanoma_mask <- getTissueType(sampleNames(ccle)) %in% c("SKIN")
skcm_braf <- virtual_ic50(ccle, "PLX4720",skcm, seed=2013)
skcm_braf_vogel <- find_drug_features (skcm_braf,skcm, with.rppa=FALSE,beta_threshold=10^-3,gene.dict="vogelstein")
skcm_braf_cosmic <- find_drug_features (skcm_braf,skcm, with.rppa=FALSE,beta_threshold=10^-3,gene.dict="cosmic")

skcm_braf_vogel$df[1:10,]

pdf("plots/tcga/skcm_plx4720.pdf",width=10,height=6,useDingbats=F)
plot_features(skcm_braf_DF, "PLX4720 (braf inib)","melanoma",text.cex=.9)
dev.off()
pdf("plots/tcga/skcm_rppa_plx4720.pdf",width=10,height=6,useDingbats=F)
plot_features(skcm_braf_rppa_DF, "PLX4720 (braf inib)","melanoma",text.cex=.9)
dev.off()

skcm_braf_2 <- virtual_ic50(ccle, "SORAFENIB",skcm, seed=2013)


##########################
## MEK in crc
carcinoma_mask <- !(getTissueType(sampleNames(ccle)) %in% c("CENTRAL_NERVOUS_SYSTEM","SKIN","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
coad_mek <- virtual_ic50(ccle[,carcinoma_mask], "PD0325901",crc, seed=2013)
coad_mek_driver <- find_drug_features (coad_mek,crc, with.rppa=FALSE,beta_threshold=10^-3,driver.genes.only=TRUE)
coad_mek_cgenes <- find_drug_features (coad_mek,crc, with.rppa=FALSE,beta_threshold=10^-3,driver.genes.only=FALSE)
coad_mek_cgenes_rppa <- find_drug_features (coad_mek,crc, with.rppa=TRUE,beta_threshold=10^-3,driver.genes.only=FALSE)

coad_mek_null <- t(replicate(100, find_drug_features (coad_mek,crc, with.rppa=FALSE,beta_threshold=10^-3,num.bootstraps=5,randomize=TRUE)$metric))
coad_mek_rppa_null <- t(replicate(100, find_drug_features (coad_mek,crc, with.rppa=TRUE,beta_threshold=10^-3,num.bootstraps=5,randomize=TRUE)$metric))

pdf("plots/tcga/crc_mek_dgenes.pdf",width=10,height=6,useDingbats=F)
plot_features(coad_mek_driver, "MEK","CRC",text.cex=.9)
dev.off()
pdf("plots/tcga/crc_mek_cgenes.pdf",width=10,height=6,useDingbats=F)
plot_features(coad_mek_cgenes, "MEK","CRC",text.cex=.9)
dev.off()
pdf("plots/tcga/crc_mek_cgenes_rppa.pdf",width=10,height=6,useDingbats=F)
plot_features(coad_mek_cgenes_rppa, "MEK","CRC",text.cex=.9)
dev.off()

#########################
# erlotinib in LUNG
lung_mask <- getTissueType(sampleNames(ccle)) %in% c("LUNG")
luad_erl <- virtual_ic50(ccle[,lung_mask], "ERLOTONIB",luad, seed=2013)
luad_erl_DF <- find_drug_features (luad_erl,luad, with.rppa=FALSE,beta_threshold=10^-3)
luad_erl_rppa_DF <- find_drug_features (luad_erl,luad, with.rppa=TRUE,beta_threshold=10^-3)

pdf("plots/tcga/luad_erl.pdf",width=10,height=6,useDingbats=F)
plot_features(luad_erl_DF, "Erlotinib","LUAD",text.cex=.9)
dev.off()
pdf("plots/tcga/luad_rppa_erl.pdf",width=10,height=6,useDingbats=F)
plot_features(luad_erl_rppa_DF, "Erlotinib","LUAD",text.cex=.9)
dev.off()

###################
# breast, pik3 pathway
breast_mask <- grepl("BREAST",getTissueType(sampleNames(sanger)))
carcinoma_mask <- !(getTissueType(sampleNames(sanger)) %in% c("CENTRAL_NERVOUS_SYSTEM","SKIN","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
par(mfrow=c(2,3))
y_hat_akt_1 <- virtual_ic50(sanger[,breast_mask],"A-443654",brca,reverseDrug=TRUE)
y_hat_akt_2 <-  virtual_ic50(sanger[,breast_mask],"AKT-inhibitor-VIII",brca,reverseDrug=TRUE)
y_hat_pik3b <- virtual_ic50(sanger[,breast_mask],"AZD6482",brca,reverseDrug=TRUE)
y_hat_pi3k_1 <- virtual_ic50(sanger[,breast_mask],"GDC0941",brca,reverseDrug=TRUE)
y_hat_pi3k_2 <- virtual_ic50(sanger[,breast_mask],"NVP-BEZ235",brca,reverseDrug=TRUE)

###############
## fgfr in breast
carcinoma_mask <- !(getTissueType(sampleNames(sanger)) %in% c("CENTRAL_NERVOUS_SYSTEM","SKIN","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
y_hat_fgfr <- virtual_ic50(sanger[,carcinoma_mask],"PD-173074",brca,reverseDrug=TRUE)
brca_fgfr_vogel <- find_drug_features (y_hat_fgfr,brca, with.rppa=FALSE,beta_threshold=10^-3,gene.dict="vogelstein")
brca_fgfr_cosmic <- find_drug_features (y_hat_fgfr,brca, with.rppa=FALSE,beta_threshold=10^-3,gene.dict="cosmic")
brca_fgfr_cbio <- find_drug_features (y_hat_fgfr,brca, with.rppa=FALSE,beta_threshold=10^-3,gene.dict="cbio")
pdf("plots/brca_fgfr_vogel.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_fgfr_vogel, "PD-173074 (FGFR)","BRCA",text.cex=.9)
dev.off()
pdf("plots/brca_fgfr_cosmic.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_fgfr_cosmic, "PD-173074 (FGFR)","BRCA",text.cex=.9)
dev.off()
pdf("plots/brca_fgfr_cbio.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_fgfr_cbio, "PD-173074 (FGFR)","BRCA",text.cex=.9)
dev.off()


####################
## MTOR in breast 
y_hat_mtor <- virtual_ic50(sanger[,carcinoma_mask],"Temsirolimus",brca,reverseDrug=TRUE)
brca_mtor_dgenes <- find_drug_features (y_hat_mtor,brca, with.rppa=FALSE,beta_threshold=10^-3,driver.genes.only=TRUE)
brca_mtor_cgenes <- find_drug_features (y_hat_mtor,brca, with.rppa=FALSE,beta_threshold=10^-3,driver.genes.only=FALSE)
brca_mtor_cgenes_rppa <- find_drug_features (y_hat_mtor,brca, with.rppa=TRUE,beta_threshold=10^-3,driver.genes.only=FALSE)
pdf("plots/tcga/brca_mtor_dgenes.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_mtor_dgenes, "Temsirolimus","BRCA",text.cex=.9)
dev.off()
pdf("plots/tcga/brca_mtor_cgenes.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_mtor_cgenes, "Temsirolimus","BRCA",text.cex=.9)
dev.off()
pdf("plots/tcga/brca_mtor_cgenes_rppa.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_mtor_cgenes_rppa, "Temsirolimus","BRCA",text.cex=.9)
dev.off()

################
## MTOR 3neg bc
y_hat_mtor <- virtual_ic50(sanger[,carcinoma_mask],"Temsirolimus",brca.3neg,reverseDrug=TRUE)
brca_mtor_cgenes <- find_drug_features (y_hat_mtor,brca.3neg, with.rppa=FALSE,beta_threshold=10^-3,driver.genes.only=FALSE,min.count=1)
brca_mtor_cgenes_rppa <- find_drug_features (y_hat_mtor,brca.3neg, with.rppa=TRUE,beta_threshold=10^-3,driver.genes.only=FALSE,min.count=1)
pdf("plots/tcga/brca3neg_mtor_cgenes.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_mtor_cgenes, "Temsirolimus","BRCA",text.cex=.9)
dev.off()
pdf("plots/tcga/brca3neg_mtor_cgenes_rppa.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_mtor_cgenes_rppa, "Temsirolimus","BRCA",text.cex=.9)
dev.off()

##################
## her2/egfr in breast
carcinoma_mask <- !(getTissueType(sampleNames(ccle)) %in% c("CENTRAL_NERVOUS_SYSTEM","SKIN","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
y_hat_her2 <- virtual_ic50(ccle[,carcinoma_mask], "LAPATINIB",brca, seed=2013)

brca_her2_dgenes <- find_drug_features (y_hat_her2,brca, with.rppa=FALSE,beta_threshold=10^-3,driver.genes.only=TRUE)
brca_her2_cgenes <- find_drug_features (y_hat_her2,brca, with.rppa=FALSE,beta_threshold=10^-3,driver.genes.only=FALSE)
brca_her2_cgenes_rppa <- find_drug_features (y_hat_her2,brca, with.rppa=TRUE,beta_threshold=10^-3,driver.genes.only=FALSE)
pdf("plots/tcga/brca_her2_dgenes.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_her2_dgenes, "Lapatanib","BRCA",text.cex=.9)
dev.off()
pdf("plots/tcga/brca_her2_cgenes.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_her2_cgenes, "Lapatanib","BRCA",text.cex=.9)
dev.off()
pdf("plots/tcga/brca_her2_cgenes_rppa.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_her2_cgenes_rppa, "Lapatanib","BRCA",text.cex=.9)
dev.off()

###################
# PIK3B in breast
brca_pik3b_DF <- find_drug_features (y_hat_pik3b,brca, with.rppa=FALSE,beta_threshold=10^-3)
brca_pik3b_rppa_DF <- find_drug_features (y_hat_pik3b,brca, with.rppa=TRUE,beta_threshold=10^-3)
pdf("plots/tcga/brca_pik3b.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_pik3b_DF[1:25,], "AZD6482 (pik3b)","BRCA",length(y_hat_mtor),text.cex=.9)
dev.off()
pdf("plots/tcga/brca_rppa_pik3b.pdf",width=10,height=6,useDingbats=F)
plot_features(brca_pik3b_rppa_DF[1:25,], "AZD6482 (pik3b)","BRCA",length(y_hat_mtor),text.cex=.9)
dev.off()

