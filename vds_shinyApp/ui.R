library(shiny)

##UI inputs and outputs

shinyUI(pageWithSidebar(

headerPanel("Virtual Drug Sensitivity Pipeline"),
sidebarPanel(
checkboxGroupInput("cancer","Select one or more cancers:",choices=cancerscheck),
actionButton("selectall","Select all cancers"),
actionButton("unselect","Select none"),
selectInput("druggroup","CCLE or Sanger drugs?",choices=c("CCLE","Sanger")),
htmlOutput("drugchoice"), #list of drugs to choose from depends on if CCLE or Sanger is selected
#selectInput("drug","Choose a drug:",choices=drugs),
conditionalPanel(condition="input.tabset1=='Features Summary'",
                 selectInput("all_subset","Display all features or only most important?",choices=c("Most important","All")),
                 selectInput("sortfeat","Sort features by:",choices=c("Importance score"="freqcounts","Effect size"="beta","P-value"="pval"))
                 ),
conditionalPanel(condition="input.tabset1=='Features Summary'&&input.all_subset=='Most important'",
                 numericInput("num_feat","Number of features to display",value=10,min=1,max=100)),

#choose a gene for by-feature results
conditionalPanel(condition="input.tabset1=='Compare Feature Across Drugs'|input.tabset1=='Compare Feature Across Cancers'",
                 selectInput("gene","Choose a gene:",choices=cclegenes)),
#plotting options for results for different drugs
#conditionalPanel(condition="input.tabset1=='Compare Feature Across Drugs'",
 #                selectInput("xaxis","Plot X variable",choices=c("Effect size"="beta","Importance score"="freqcounts","P-value"="pval")),
  #               selectInput("yaxis","Plot Y variable",choices=c("Effect size"="beta","Importance score"="freqcounts","P-value"="pval")),
   #              actionButton("plotdrugs","Plot Drugs")),
conditionalPanel(condition="input.tabset1=='Compare Feature Across Drugs'",
                 selectInput("drugsort","Sort table by:",choices=c("Importance score"="freqcounts","Effect size"="beta","P-Value"="pval")),
                 numericInput("num_drugs","Number of drugs to show in table (max 13 for paper results, 154 for all)",value=10,min=1,max=154)
                 ),
                 checkboxInput("paperonly","Show only most recent results included in paper"),
#This was never actually implemented - show results for drugs with same target
#conditionalPanel(condition="input.tabset1=='Compare Feature Across Cancers'",
#                 actionButton("reldrugs","See results for this gene for other related drugs")),

##KEGG pathway visualization
conditionalPanel(condition="input.tabset1=='KEGG Pathway View'",
                 radioButtons("keggmethod","Choose from main KEGG cancer pathways or view a list of the most enriched pathways for this cancer and drug:",choices=c("Main Pathways"="main","Most Enriched"="calc")),
                 conditionalPanel(condition="input.keggmethod=='main'",
                                  selectInput("keggid","Choose KEGG pathway from main cancer pathways to display",choices=keggpathways)
                                  ),
                 conditionalPanel(condition="input.keggmethod=='calc'",
                                uiOutput("topmapchoices")),
                 selectInput("highlight","Highlight nodes based on:",choices=c("Model effect"="beta","Importance score"="freqcounts","Frequency in TCGA population"="freqevents"))
                 #actionButton("makepath","Make and show Pathview")
                 )),
mainPanel(
tabsetPanel(
tabPanel( "Model performance",
h5("Performance metrics"),
#textOutput("auccelllines"),
conditionalPanel(condition="input.drug.length>0",
                 tableOutput("perf"),
                  uiOutput("perfall")),
uiOutput("allperfplot"),
uiOutput("cancerperfplots")),
tabPanel("Features Summary",
         uiOutput("featuresummary")),
tabPanel("Compare Feature Across Cancers",
         uiOutput("drugname"),
         uiOutput("feat_tables")
         ),
tabPanel("Compare Feature Across Drugs",
         #conditionalPanel(condition="input.cancer.length>0",
          #                plotOutput("featDrugPlots")),
         conditionalPanel(condition= "input.cancer.indexOf('blca')!=-1&&input.paperonly!=1",
                          h5("Bladder urothelial carcinoma"),
                          uiOutput("blcatables")),
         conditionalPanel(condition="input.cancer.indexOf('kirc')!=-1",
                          h5("Kidney renal cell carcinoma"),
                          uiOutput("kirctables")
                          ),
         conditionalPanel(condition="input.cancer.indexOf('gbm')!=-1",
                          h5("Glioblastoma multiforme"),
                          uiOutput("gbmtables")),
         conditionalPanel(condition="input.cancer.indexOf('laml')!=-1&&input.paperonly!=1",
                          h5("Acute myeloid leukemia"),
                          uiOutput("lamltables")),
         conditionalPanel(condition="input.cancer.indexOf('crc')!=-1",
                          h5("Colorectal cancer"),
                          uiOutput("crctables")),
         conditionalPanel(condition="input.cancer.indexOf('ov')!=-1",
                          h5("Ovarian serous cystadenocarcinoma"),
                          uiOutput("ovtables")),
         conditionalPanel(condition= "input.cancer.indexOf('ucec')!=-1",
                          h5("Uterine corpus endometrial carcinoma"),
                          uiOutput("ucectables")),
         conditionalPanel(condition="input.cancer.indexOf('prad')!=-1&&input.paperonly!=1",
                          h5("Prostate adenocarcinoma"),
                          uiOutput("pradtables")),
         conditionalPanel(condition="input.cancer.indexOf('skcm')!=-1",
                          h5("Skin cutaneous melanoma"),
                          uiOutput("skcmtables")),
         conditionalPanel(condition="input.cancer.indexOf('lusc')!=-1",
                          h5("Lung squamous cell carcinoma"),
                          uiOutput("lusctables")),
         conditionalPanel(condition="input.cancer.indexOf('luad')!=-1",
                          h5("Lung adenocarcinoma"),
                          uiOutput("luadtables")),
         conditionalPanel(condition="input.cancer.indexOf('brca')!=-1",
                          h5("Breast invasive carcinoma"),
                          uiOutput("brcatables")),
         conditionalPanel(condition="input.cancer.indexOf('brca.3neg')!=-1",
                          h5("Triple negative breast cancer"),
                          uiOutput("brca.3negtables")),
         conditionalPanel(condition="input.cancer.indexOf('brca.erpr')!=-1",
                          h5("ER/PR-positive breast cancer"),
                          uiOutput("brca.erprtables")),
         conditionalPanel(condition="input.cancer.indexOf('brca.her2')!=-1",
                          h5("HER2-positive breast cancer"),
                          uiOutput("brca.her2tables")),
         conditionalPanel(condition="input.cancer.indexOf('stad')!=-1",
                          h5("Stomach adenocarcinoma"),
                          uiOutput("stadtables"))
         ),
tabPanel("KEGG Pathway View",
         helpText("Only the first cancer selected is used to generate the pathway maps. The image file is saved as keggID.cancer_Drug.png"),
          uiOutput("pathwaymap")
         ),
id="tabset1")
)
))