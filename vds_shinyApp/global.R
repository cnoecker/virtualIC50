library(pathview)
library(org.Hs.eg.db)
library(KEGGREST)
options("scipen"=100,digits=6) ##avoid scientific notation


###### Define drugs and cancers

drugs=c("ccle_AEW541","sanger_Etoposide","ccle_CRIZOTINIB","ccle_DEMETHOXYGELDANAMYCIN","sanger_FH535",
        "sanger_FTI-277","ccle_ERLOTONIB","sanger_GDC-0449","ccle_IRINOTECAN","sanger_GDC0941",
        "ccle_L685458","sanger_Gefitinib","ccle_LAPATINIB","sanger_Gemcitabine","ccle_LBW242",
        "sanger_GNF-2","ccle_NILOTINIB","sanger_GSK269962A","ccle_NUTLIN3","sanger_GSK-650394",
        "ccle_PACLITAXEL","sanger_GW-441756","ccle_PANOBINOSTAT","sanger_GW843682X","ccle_PD0325901",
        "sanger_Imatinib","ccle_PD0332991","sanger_IPA-3","ccle_PHA665752","sanger_JNK-9L","ccle_PLX4720",
        "sanger_JNK-Inhibitor-VIII","ccle_RAF265","sanger_JW-7-52-1","ccle_SARACATINIB","sanger_KIN001-135",
        "ccle_SELUMETINIB",            "sanger_KU-55933",
        "ccle_SORAFENIB",              "sanger_Lapatinib",
        "ccle_TAE684",                 "sanger_Lenalidomide",
        "ccle_TKI258",                 "sanger_LFM-A13",
        "ccle_TOPOTECAN",              "sanger_Metformin",
        "ccle_VANDETANIB",             "sanger_Methotrexate",
        "sanger_A-443654",             "sanger_MG-132",
        "sanger_A-769662",             "sanger_Midostaurin",
        "sanger_A-770041",             "sanger_Mitomycin-C",
        "sanger_ABT-263",              "sanger_MK-2206",
        "sanger_ABT-888",              "sanger_MS-275",
        "sanger_AICAR",                "sanger_Nilotinib",
        "sanger_AKT-inhibitor-VIII",   "sanger_NSC-87877",
        "sanger_AMG-706",              "sanger_NU-7441",
        "sanger_AP-24534",             "sanger_Nutlin-3a",
        "sanger_AS601245",             "sanger_NVP-BEZ235",
        "sanger_ATRA",                 "sanger_NVP-TAE684",
        "sanger_AUY922",               "sanger_Obatoclax-Mesylate",
        "sanger_Axitinib",             "sanger_OSI-906",
        "sanger_AZ628",                "sanger_PAC-1",
        "sanger_AZD-0530",             "sanger_Paclitaxel",
        "sanger_AZD-2281",             "sanger_Parthenolide",
        "sanger_AZD6244",              "sanger_Pazopanib",
        "sanger_AZD6482",              "sanger_PD-0325901",
        "sanger_AZD7762",              "sanger_PD-0332991",
        "sanger_AZD8055",              "sanger_PD-173074",
        "sanger_BAY-61-3606",          "sanger_PF-02341066",
        "sanger_Bexarotene",           "sanger_PF-562271",
        "sanger_BI-2536",              "sanger_PHA-665752",
        "sanger_BIBW2992",             "sanger_PLX-4720",
        "sanger_Bicalutamide",         "sanger_Pyrimethamine",
        "sanger_BI-D1870",             "sanger_QS11",
        "sanger_BIRB-0796",            "sanger_Rapamycin",
        "sanger_Bleomycin",            "sanger_RDEA119",
        "sanger_BMS-509744",           "sanger_RO-3306",
        "sanger_BMS-536924",           "sanger_Roscovitine",
        "sanger_BMS-754807",           "sanger_Salubrinal",
        "sanger_Bortezomib",           "sanger_SB-216763",
        "sanger_Bosutinib",            "sanger_SB590885",
        "sanger_Bryostatin-1",         "sanger_Shikonin",
        "sanger_BX-795",               "sanger_SL-0101-1",
        "sanger_Camptothecin",         "sanger_Sorafenib",
        "sanger_CEP-701",              "sanger_S-Trityl-L-cysteine",
        "sanger_CGP-082996",           "sanger_Sunitinib",
        "sanger_CGP-60474",            "sanger_Temsirolimus",
        "sanger_CHIR-99021",           "sanger_Thapsigargin",
        "sanger_CI-1040",              "sanger_Tipifarnib",
        "sanger_Cisplatin",            "sanger_Vinblastine",
        "sanger_CMK",                  "sanger_Vinorelbine",
        "sanger_Cyclopamine",          "sanger_Vorinostat",
        "sanger_Cytarabine",           "sanger_VX-680",
        "sanger_Dasatinib",            "sanger_VX-702",
        "sanger_DMOG",                 "sanger_WH-4-023",
        "sanger_Docetaxel",            "sanger_WZ-1-84",
        "sanger_Doxorubicin",          "sanger_X17-AAG",
        "sanger_Elesclomol",           "sanger_X681640",
        "sanger_Embelin",              "sanger_XMD8-85",
        "sanger_Epothilone-B",         "sanger_Z-LLNle-CHO",
        "sanger_Erlotinib",            "sanger_ZM-447439"
)
drugs=sort(drugs)
drugs2=gsub("sanger_","",drugs)
drugs2=gsub("ccle_","",drugs2)


ccledrugs=drugs[grepl("ccle",drugs)]
ccledrugs=gsub("ccle_","",ccledrugs)

sangerdrugs=drugs[grepl("sanger",drugs)]
sangerdrugs=gsub("sanger_","",sangerdrugs)

paperdrugs=c("AZD-2281","Bortezomib","Cisplatin","Doxorubicin","ERLOTONIB","LAPATINIB","PACLITAXEL","PLX4720","SELUMETINIB",
             "SORAFENIB","Sunitinib","Temsirolimus","TKI258")
papercancers=c("kirc","gbm","crc","ov","ucec","skcm","lusc","luad","brca","brca.3neg","brca.erpr","brca.her2","stad")

cancers=c("blca","kirc","gbm","laml","crc","ov","ucec","prad","skcm","lusc","luad","brca")
allcancers=c(cancers,papercancers[10:13])

cancerscheck=c("Bladder urothelial carcinoma (BLCA)"="blca",
               "Kidney renal clear cell carcinoma (KIRC)"="kirc",
               "Glioblastoma multiforme (GBM)"="gbm",
               "Acute myeloid leukemia (LAML)"="laml",
               "Colorectal cancer (CRC)"="crc",
               "Ovarian serous cystadenocarcinoma (OV)"="ov",
               "Uterine corpus endometrial carcinoma (UCEC)"="ucec",
               "Prostate adenocarcinoma (PRAD)"="prad",
               "Skin cutaneous melanoma (SKCM)"="skcm",
               "Lung squamous cell carcinoma (LUSC)"="lusc",
               "Lung adenocarcinoma (LUAD)"="luad",
               "Breast invasive carcinoma (BRCA)"="brca",
               "Triple negative breast cancer (only for select drugs)"="brca.3neg",
               "ER/PR-positive breast cancer (only for select drugs)"="brca.erpr",
               "HER2-positve breast cancer (only for select drugs)"="brca.her2",
               "Stomach adenocarcinoma (STAD) (only for select drugs)"="stad")
cancerlabels=c("Bladder urothelial carcinoma (BLCA)",
               "Kidney renal clear cell carcinoma (KIRC)",
               "Glioblastoma multiforme (GBM)",
               "Acute myeloid leukemia (LAML)",
               "Colorectal cancer (CRC)",
               "Ovarian serous cystadenocarcinoma (OV)",
               "Uterine corpus endometrial carcinoma (UCEC)",
               "Prostate adenocarcinoma (PRAD)",
               "Skin cutaneous melanoma (SKCM)",
               "Lung squamous cell carcinoma (LUSC)",
               "Lung adenocarcinoma (LUAD)",
               "Breast invasive carcinoma (BRCA)",
               "Triple negative breast cancer (only for select drugs)",
               "ER/PR-positive breast cancer (only for select drugs)",
               "HER2-positve breast cancer (only for select drugs)",
               "Stomach adenocarcinoma (STAD) (only for select drugs)")

keggpathways=c("Cell cycle"="04110","MAPK signaling pathway"="04010", "p53 signaling pathway"="04115","Apoptosis"="04210","Pathways in cancer"="05200",
  "Transcriptional misregulation in cancer"="05202","MicroRNAs in cancer"="05206",
  "Proteoglycans in cancer"="05205","Chemical carcinogenesis"="05204","Viral carcinogenesis"="05203",
  "Colorectal cancer"="05210","Endometrial cancer"="05213","Glioma"="05214","Acute myeloid leukemia"="05221",
  "Melanoma"="05218","Renal cell carcinoma"="05211","Bladder cancer"="05219","Prostate cancer"="05215",
  "Small cell lung cancer"="05222")

### Model results - see resultsRead.R
load("/gluster/home/cnoecker/RDataVdsApp/vdsResults.rda")

######KEGG info - keggPaths_FisherTests.R
load("/gluster/home/cnoecker/RDataVdsApp/keggpathways.rda") 
load("/gluster/home/cnoecker/RDataVdsApp/fisherpathresults.rda")


##Genes and genetic features

cclefeatures=list(NULL)
for(i in 1:length(cancers)){
  cclefeatures[[i]]=sort(results[[1]][[i]]$df$genes)
}

cclefeatures2=sort(unique(unlist(cclefeatures)))

sangerfeatures=list(NULL)
for(i in 1:length(cancers)){
  sangerfeatures[[i]]=sort(results[[25]][[i]]$df$genes)
}
sangerfeatures2=sort(unique(unlist(sangerfeatures)))

paperfeatures=list(NULL)
for(i in 1:length(allcancers)){
  paperfeatures[[i]]=sort(results[[40]][[i]]$df$genes)
}
paperfeatures2=sort(unique(unlist(paperfeatures)))

cclegenes=c()
for(i in 1:length(cclefeatures)){
  for(j in 1:length(cclefeatures[[i]])){
    genevec=c()
    genevec=strsplit(as.character(cclefeatures[[i]][j]),"::")
    cclegenes=c(cclegenes,genevec)
    }
  }
cclegenes=sort(unique(unlist(cclegenes)))

cclegenes=gsub("mutORamp","",cclegenes)
cclegenes=gsub("delORmut","",cclegenes)
cclegenes=gsub("_","",cclegenes)
cclegenes=gsub("amp","",cclegenes)
cclegenes=gsub("mut","",cclegenes)
cclegenes=gsub("del","",cclegenes)
cclegenes=unique(cclegenes)

sangergenes=c()
for(i in 1:length(sangerfeatures)){
  for(j in 1:length(sangerfeatures[[i]])){
    genevec=c()
    genevec=strsplit(as.character(sangerfeatures[[i]][j]),"::")
    sangergenes=c(sangergenes,genevec)
  }
}
sangergenes=sort(unique(unlist(sangergenes)))

sangergenes=gsub("mutORamp","",sangergenes)
sangergenes=gsub("delORmut","",sangergenes)
sangergenes=gsub("_","",sangergenes)
sangergenes=gsub("amp","",sangergenes)
sangergenes=gsub("mut","",sangergenes)
sangergenes=gsub("del","",sangergenes)
sangergenes=unique(sangergenes)

cclefeats=c()
for(i in 1:length(cclegenes)){
  grep1=grep(paste("^",cclegenes[i],"_",sep=""),cclefeatures2)
  grep2=grep(paste(":",cclegenes[i],"_",sep=""),cclefeatures2)
  cclefeats[i]=length(c(grep1,grep2))
}
maxfeats=max(cclefeats) #### 11 different features for 1 gene

# sangerfeats=c()
# for(i in 1:length(sangergenes)){
#   grep1=grep(paste("^",sangergenes[i],"_",sep=""),sangerfeatures2)
#   grep2=grep(paste(":",sangergenes[i],"_",sep=""),sangerfeatures2)
#   sangerfeats[i]=length(c(grep1,grep2))
# }
# maxfeats2=max(sangerfeats)




############## Show results for drugs with the same target - never fully implemented

# targetinfo=read.csv("/gluster/home/cnoecker/drugtargets.csv")
# head(targetinfo)
# table(targetinfo$Target)
# drugs2=gsub("ccle_","",drugs)
# drugs2=gsub("sanger_","",drugs2)
# 
# targets=rep(NA,length(drugs2))
# for(i in 1:length(drugs2)){
#   ind=match(drugs2[i],targetinfo$Drug)
#     targets[i]=as.character(targetinfo$Target[ind])
# }
# 
# write.csv(cbind(drugs2,targets),file="/gluster/home/cnoecker/targets.csv")


# targetinfo=read.csv("/gluster/home/cnoecker/targets.csv")
# targets=targetinfo$targets
# targets=gsub(",","",targets)

