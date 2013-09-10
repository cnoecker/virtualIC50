library(shiny)
library(RCircos)
library(reshape)

shinyServer(function(input, output,clientData,session) {
  
  observe({
      if(input$selectall)
      updateCheckboxGroupInput(session,inputId="cancer",choices=cancerscheck,selected=cancerlabels)
    
      if(input$unselect)
      updateCheckboxGroupInput(session,inputId="cancer",choices=cancerscheck)
  })
  
  DrugCancer=function(x){
    results[[drugi()]][[x]]
  }
  
  ##Cancer and drug indices
  indices=reactive({
    sapply(1:numcanc(),function(x){match(input$cancer[x],allcancers)})
  })
  
  numcanc=reactive({
    length(input$cancer)
  })
  
  drugi=reactive({
    if(!is.null(input$drug)){
      if(identical(input$druggroup,"CCLE"))
        grep(input$drug,ccledrugs) else
          if(identical(input$druggroup,"Sanger"))
            grep(input$drug,sangerdrugs)+24}
  })
  
  #Reactive drug choices
  output$drugchoice=renderUI({
    choices=c()
    if(input$druggroup=="CCLE"){
      if(any(input$cancer %in% papercancers[10:13]))
        choices=paperdrugs[paperdrugs %in% ccledrugs]
      else choices=ccledrugs
    } else if(input$druggroup=="Sanger"){
      if(any(input$cancer %in% papercancers[10:13]))
        choices=paperdrugs[paperdrugs %in% sangerdrugs]
      else choices=sangerdrugs
    }
    selectInput("drug","Choose a drug",choices=choices)
  })
  
  ##Feature lists to use depending on subset of results
  features=function(drugi){
    if(drugs2[drugi] %in% paperdrugs) return(list(paperfeatures)) else if(identical(input$druggroup,"Sanger"))
      return(list(sangerfeatures)) else if(identical(input$druggroup,"CCLE")) return(list(cclefeatures))
  }
  
  features2=function(drugi){
    if(drugs2[drugi][1] %in% paperdrugs) return(list(paperfeatures2)) else if(identical(input$druggroup,"Sanger"))
      return(list(sangerfeatures2)) else if(identical(input$druggroup,"CCLE")) return(list(cclefeatures2))
  }
  
  #Selected features
  genefeati=reactive({
    c(grep(paste("^",input$gene,"_",sep=""),paperfeatures2),grep(paste(":",input$gene,"_",sep=""),
                                                                                           paperfeatures2))
  })
  
  myfeat=reactive({
    paperfeatures2[genefeati()]
  })
  
  ##filename prefix for loading image files
  filedate=reactive({
    if(input$druggroup=="CCLE") "all_2013_08_03/ccle_" else if(input$druggroup=="Sanger") "all_2013_06_21/sanger_"
  })
  
  ##Performance metrics 
#   output$auccelllines=renderText({
#     if(!is.null(input$drug))
#       paste("AUC for model of ",input$drug,"over all cell lines:",round(auc_allcelllines[drugi()],digits=4))
#   })
   
  output$perf=renderTable({
    if(numcanc()>0&!is.null(input$drug)){
      perftable=data.frame(Cancer=rep(NA,numcanc()),auc=rep(NA,numcanc()),rho=rep(NA,numcanc()),fdr=rep(NA,numcanc()),N=rep(NA,numcanc())) 
      for(i in 1:numcanc()){
        fdr<- if(input$druggroup=="CCLE"&is.na(match(input$drug,paperdrugs))){DrugCancer(indices()[i])$fdr}else {NA}
        perftable[i,]=c(input$cancer[i],round(DrugCancer(indices()[i])$metric[2],digits=4),
                        round(DrugCancer(indices()[i])$metric[1],digits=4),round(fdr,digits=4),DrugCancer(indices()[i])$N)
      }
      return(perftable)
    }
  },include.rownames=F,digits=c(0,4,4,4,4,0),align=rep("c",6))

  output$perfall=renderUI({
    if(!is.null(input$drug)){
    if(input$drug %in% paperdrugs){
    list(h5("Performance metrics across all cell lines"),
    tableOutput("perfalltab"))
}
    }
  })
  
  output$perfalltab=renderTable({
    paperperftable[paperperftable$drug==input$drug,]
  },include.rownames=F,digits=c(0,3,3,3,3,3),align=rep("c",6))
  
  output$allperfplot=renderUI({
    if(!is.null(input$drug))
    if(input$drug %in% paperdrugs)
    imageOutput("paperperf")
  })
  
  output$paperperf=renderImage({
    filename=normalizePath(file.path(paste('/gluster/home/jguinney/projects/virtualIC50/results/paper_2013_08_30/perfPlots/',input$drug,'.pdf',sep='')))
    list(src=filename,width=500)
  },deleteFile=FALSE)
  
  ##Performance plots
  imagefile=function(canceri){
    filename=normalizePath(file.path(paste('/gluster/home/jguinney/projects/virtualIC50/results/',filedate(),
                                           input$drug,'/',allcancers[canceri],'_perf.pdf',sep='')))
  return(filename)
  }
  
  ####FOR MONDAY - what is happening?????
  papercanc=reactive({
    if(!is.null(input$drug)&length(input$cancer)>0){
      if(input$drug %in% paperdrugs){
        papercanc=input$cancer[!input$cancer %in% papercancers]
      } else papercanc=input$cancer
    }
  })
  
  output$cancerperfplots=renderUI({
    outputlistperf=list(NULL)
    if(!is.null(input$drug)&length(papercanc())>0){
    for(i in 1:length(papercanc())){
      outputlistperf[[2*i-1]]<-{
        cancername=paste("perf",i,sep="")
        textOutput(cancername)
          }
      outputlistperf[[2*i]]<-{
        plotname=paste("pplot",i,sep="")
        imageOutput(plotname,height="auto") 
      }
    }
    do.call(tagList,outputlistperf)
    }
  })

  for(k in 1:length(allcancers)){
    local({
      my_k=k
      cancername=paste("perf",my_k,sep="")
      output[[cancername]]=renderText({
        cancerlabels[match(papercanc()[my_k],allcancers)]
      })
      plotname=paste("pplot",my_k,sep="")
      output[[plotname]]<-renderImage({
        list(src=imagefile(match(papercanc()[my_k],allcancers)),width=500)
      },deleteFile=FALSE) 
      
      plotnamef=paste("fplot",my_k,sep="")
      output[[plotnamef]]=renderImage({
        canceri=match(input$cancer[my_k],allcancers)
        list(src=fplotfile(canceri),width=700)
      },deleteFile=F)
      cancernamef=paste("feat",my_k,sep="")
      output[[cancernamef]]=renderTable({
        canceri=match(input$cancer[my_k],allcancers)
        if((identical("Most important",input$all_subset)) & isTRUE(input$num_feat<length(features(match(input$drug,drugs2))[[1]][[canceri]])))
          effects_tab(canceri)[1:input$num_feat,] else if (identical("All",input$all_subset)|input$num_feat>length(features(match(input$drug,drugs2))[[1]][[canceri]]))
            effects_tab(canceri)
      },digits=c(0,0,2,4,3,3,4),align=rep("c",7),include.rownames=F)
    })
  }

  ##Load feature (bubble) plots
  fplotfile=function(canceri){
    if(!is.na(match(input$drug,paperdrugs))&!is.na(match(allcancers[canceri],papercancers))){
      filename=normalizePath(file.path(paste('/gluster/home/jguinney/projects/virtualIC50/results/paper_2013_08_30/bubblePlots/',input$drug,'_',allcancers[canceri],'.pdf',sep='')))
}else filename=normalizePath(file.path(paste('/gluster/home/jguinney/projects/virtualIC50/results/',filedate(),
                                           input$drug,'/',allcancers[canceri],'.pdf',sep='')))
    return(filename)
  }
  
  output$featuresummary=renderUI({
    outputlistfeat=list(NULL)
    if(!is.null(input$drug)&length(input$cancer)>0){
      for(i in 1:length(input$cancer)){
        outputlistfeat[[2*i-1]]<-{
          plotnamef=paste("fplot",i,sep="")
          imageOutput(plotnamef)
        }
        outputlistfeat[[2*i]]<-{
          cancernamef=paste("feat",i,sep="")
          tableOutput(cancernamef)
        }
      }
      do.call(tagList,outputlistfeat)
    }
  })
    
  ##Effects table function for 
  effects_tab=function(x){
    if(identical(input$druggroup,"CCLE")){
      effectsdf=data.frame(Feature=DrugCancer(x)$df[,1],posFreq=DrugCancer(x)$df[,3],effect=DrugCancer(x)$df[,4],CountsScore=DrugCancer(x)$df[,5],freqFeature=DrugCancer(x)$df[,7],pVal=DrugCancer(x)$df[,8])
    }
    else if(identical(input$druggroup,"Sanger")){
      if(input$drug %in% paperdrugs & allcancers[x] %in% papercancers){
        effectsdf=data.frame(Feature=DrugCancer(x)$df[,1],posFreq=DrugCancer(x)$df[,3],effect=DrugCancer(x)$df[,4],CountsScore=DrugCancer(x)$df[,5],freqFeature=DrugCancer(x)$df[,7],pVal=DrugCancer(x)$df[,8])
      } else{
      neweffect=c()
      for(i in 1:length(DrugCancer(x)$df$genes)){
        neweffect[i]=subset(DrugCancer(x)$effect,!is.na(match(names(DrugCancer(x)$effect),DrugCancer(x)$df$genes[i],nomatch=NA)))
      }
      effectsdf=data.frame(Feature=DrugCancer(x)$df[,1],
                           posFreq=DrugCancer(x)$df[,3],
                           effect=neweffect,
                           CountsScore=DrugCancer(x)$df[,4],
                           freqFeature=DrugCancer(x)$df[,6],
                           pVal=DrugCancer(x)$df[,7])
      }
    }
    if(identical(input$sortfeat,"freqcounts")) return(effectsdf[order(effectsdf[,4],decreasing=T),]) else
      if(identical(input$sortfeat,"beta")) return(effectsdf[order(abs(effectsdf[,3]),decreasing=T),]) else
        if(identical(input$sortfeat,"pval")) return(effectsdf[order(effectsdf[,6],decreasing=F),])
  }
  
  
output$drugname=renderUI({
    h5(paste("Results for ",input$drug,"and",input$gene,"across all cancers\n",sep=" "))
  })
  
##rendering of Feature Details for each cancer  (see end of script for output contents)
  output$feat_tables=renderUI({
    numgene=length(genefeati())
    outputlist=list(NULL)
    for(i in 1:numgene){
      outputlist[[2*i-1]]<-{
        textname=paste("header",i,sep="")
        textOutput(textname)
      }
      outputlist[[2*i]]<-{
        tablename=paste("table",i,sep="")
        tableOutput(tablename) 
      }
    }
    do.call(tagList,outputlist)
  })
    
##Feature Details for top drugs
  output$blcatables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[1]])
    if(any(testgene)){
      outputlist2=list(NULL)
      numgene=length(genefeati())
      for(i in 1:numgene){
        outputlist2[[2*i-1]]<-{
          blcatextname=paste("blcafreq",i,sep="")
          textOutput(blcatextname)
        }
        outputlist2[[2*i]]<-{
          blcatablename=paste("blcatable",i,sep="")
          tableOutput(blcatablename)
        }
      }
      do.call(tagList,outputlist2)
    } else{
      textOutput("noresults")
    }  
  })
  output[["noresults"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$kirctables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[2]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist2a=list(NULL)
      for(i in 1:numgene){
        outputlist2a[[2*i-1]]<-{
          kirctextname=paste("kircfreq",i,sep="")
          textOutput(kirctextname)
        }
        outputlist2a[[2*i]]<-{
          kirctablename=paste("kirctable",i,sep="")
          tableOutput(kirctablename)
        }
      }
      do.call(tagList,outputlist2a)
    } else{
      textOutput("noresults2")
    }
  })
  output[["noresults2"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$gbmtables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[3]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist3=list(NULL)
      for(i in 1:numgene){
        outputlist3[[2*i-1]]<-{
          gbmtextname=paste("gbmfreq",i,sep="")
          textOutput(gbmtextname)
        }
        outputlist3[[2*i]]<-{
          gbmtablename=paste("gbmtable",i,sep="")
          tableOutput(gbmtablename)
        }
      }
      do.call(tagList,outputlist3)
    } else
      textOutput("noresults3")
  })
  output[["noresults3"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$lamltables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[4]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist4=list(NULL)
      for(i in 1:numgene){
        outputlist4[[2*i-1]]<-{
          lamltextname=paste("lamlfreq",i,sep="")
          textOutput(lamltextname)
        }
        outputlist4[[2*i]]<-{
          lamltablename=paste("lamltable",i,sep="")
          tableOutput(lamltablename)
        }
      }
      do.call(tagList,outputlist4)
    } else
      textOutput("noresults4")
  })
  output[["noresults4"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$crctables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[5]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist5=list(NULL)
      for(i in 1:numgene){
        outputlist5[[2*i-1]]<-{
          crctextname=paste("crcfreq",i,sep="")
          textOutput(crctextname)
        }
        outputlist5[[2*i]]<-{
          crctablename=paste("crctable",i,sep="")
          tableOutput(crctablename)
        }
      }
      do.call(tagList,outputlist5)
    } else textOutput("noresults5")
  })
  output[["noresults5"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$ovtables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[6]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist6=list(NULL)
      for(i in 1:numgene){
        outputlist6[[2*i-1]]<-{
          ovtextname=paste("ovfreq",i,sep="")
          textOutput(ovtextname)
        }
        outputlist6[[2*i]]<-{
          ovtablename=paste("ovtable",i,sep="")
          tableOutput(ovtablename)
        }
      }
      do.call(tagList,outputlist6)
    } else textOutput("noresults6")
  })
  output[["noresults6"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$ucectables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[7]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist7=list(NULL)
      for(i in 1:numgene){
        outputlist7[[2*i-1]]<-{
          ucectextname=paste("ucecfreq",i,sep="")
          textOutput(ucectextname)
        }
        outputlist7[[2*i]]<-{
          ucectablename=paste("ucectable",i,sep="")
          tableOutput(ucectablename)
        }
      }
      do.call(tagList,outputlist7)
    } else textOutput("noresults7")
  })
  output[["noresults7"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$pradtables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[8]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist8=list(NULL)
      for(i in 1:numgene){
        outputlist8[[2*i-1]]<-{
          pradtextname=paste("pradfreq",i,sep="")
          textOutput(pradtextname)
        }
        outputlist8[[2*i]]<-{
          pradtablename=paste("pradtable",i,sep="")
          tableOutput(pradtablename)
        }
      }
      do.call(tagList,outputlist8)
    } else textOutput("noresults8")
  })
  output[["noresults8"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$skcmtables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[9]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist9=list(NULL)
      for(i in 1:numgene){
        outputlist9[[2*i-1]]<-{
          skcmtextname=paste("skcmfreq",i,sep="")
          textOutput(skcmtextname)
        }
        outputlist9[[2*i]]<-{
          skcmtablename=paste("skcmtable",i,sep="")
          tableOutput(skcmtablename)
        }
      }
      do.call(tagList,outputlist9)
    } else textOutput("noresults9")
  })
  output[["noresults9"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$lusctables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[10]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist10=list(NULL)
      for(i in 1:numgene){
        outputlist10[[2*i-1]]<-{
          lusctextname=paste("luscfreq",i,sep="")
          textOutput(lusctextname)
        }
        outputlist10[[2*i]]<-{
          lusctablename=paste("lusctable",i,sep="")
          tableOutput(lusctablename)
        }
      }
      do.call(tagList,outputlist10)
    } else textOutput("noresults10")
  })
  output[["noresults10"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$luadtables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[11]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist11=list(NULL)
      for(i in 1:numgene){
        outputlist11[[2*i-1]]<-{
          luadtextname=paste("luadfreq",i,sep="")
          textOutput(luadtextname)
        }
        outputlist11[[2*i]]<-{
          luadtablename=paste("luadtable",i,sep="")
          tableOutput(luadtablename)
        }
      }
      do.call(tagList,outputlist11)
    } else textOutput("noresults11")
  })
  output[["noresults11"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$brcatables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[12]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist12=list(NULL)
      for(i in 1:numgene){
        outputlist12[[2*i-1]]<-{
          brcatextname=paste("brcafreq",i,sep="")
          textOutput(brcatextname)
        }
        outputlist12[[2*i]]<-{
          brcatablename=paste("brcatable",i,sep="")
          tableOutput(brcatablename)
        }
      }
      do.call(tagList,outputlist12)
    } else textOutput("noresults12")
  })
  output[["noresults12"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$brca.3negtables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[13]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist13=list(NULL)
      for(i in 1:numgene){
        outputlist13[[2*i-1]]<-{
          brca.3negtextname=paste("brca.3negfreq",i,sep="")
          textOutput(brca.3negtextname)
        }
        outputlist13[[2*i]]<-{
          brca.3negtablename=paste("brca.3negtable",i,sep="")
          tableOutput(brca.3negtablename)
        }
      }
      do.call(tagList,outputlist13)
    } else textOutput("noresults13")
  })
  output[["noresults13"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$brca.erprtables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[14]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist14=list(NULL)
      for(i in 1:numgene){
        outputlist14[[2*i-1]]<-{
          brca.erprtextname=paste("brca.erprfreq",i,sep="")
          textOutput(brca.erprtextname)
        }
        outputlist14[[2*i]]<-{
          brca.erprtablename=paste("brca.erprtable",i,sep="")
          tableOutput(brca.erprtablename)
        }
      }
      do.call(tagList,outputlist14)
    } else textOutput("noresults14")
  })
  output[["noresults14"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$brca.her2tables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[15]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist15=list(NULL)
      for(i in 1:numgene){
        outputlist15[[2*i-1]]<-{
          brca.her2textname=paste("brca.her2freq",i,sep="")
          textOutput(brca.her2textname)
        }
        outputlist15[[2*i]]<-{
          brca.her2tablename=paste("brca.her2table",i,sep="")
          tableOutput(brca.her2tablename)
        }
      }
      do.call(tagList,outputlist15)
    } else textOutput("noresults15")
  })
  output[["noresults15"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  output$stadtables=renderUI({
    testgene=grepl(input$gene,paperfeatures[[16]])
    if(any(testgene)){
      numgene=length(genefeati())
      outputlist16=list(NULL)
      for(i in 1:numgene){
        outputlist16[[2*i-1]]<-{
          stadtextname=paste("stadfreq",i,sep="")
          textOutput(stadtextname)
        }
        outputlist16[[2*i]]<-{
          stadtablename=paste("stadtable",i,sep="")
          tableOutput(stadtablename)
        }
      }
      do.call(tagList,outputlist16)
    } else textOutput("noresults16")
  })
  output[["noresults16"]]<-renderText({
    c("No results for this cancer and gene")
  })
  
  
  
##KEGG Pathway output
  
  #most enriched pathways
  topmaps=reactive({
    canceri=match(input$cancer[1],allcancers)
    ortable=fisherresults[[drugi()]][[canceri]]
    if(length(ortable[,2][as.numeric(ortable[,2])>2])>15) return(ortable[as.numeric(ortable[,2])>2,]) else return(ortable[1:15,])
  })
  
  #select from most enriched pathways
  output$topmapchoices=renderUI({
    if(input$keggmethod=="calc"){
      paths2=c()
      for(i in 1:length(topmaps()[,1])){
        paths2[i]=keggpaths2[keggpaths2==topmaps()[i,1]]
        names(paths2)[i]=names(keggpaths2[keggpaths2==topmaps()[i,1]])
      }
      selectInput("keggid2","Choose from the most enriched pathways:",choices=paths2)
    }
  })
  
  #name of pathway selected
  pathwayname=reactive({
    idinput=ifelse(input$keggmethod=="main",input$keggid,input$keggid2)
    pathwayname=names(keggpaths2[keggpaths2==idinput])
  })
  
  #generate list of values associated with gene IDs for pathview function
  pathwaygenedata=reactive({
    gene=c()
    value=c()
    canceri=match(input$cancer[1],allcancers)
    for(i in 1:length(features(drugi())[[1]][[canceri]])){
      genevec=c()
      genevec=strsplit(as.character(features(drugi())[[1]][[canceri]][i]),"::")
      genevec=unlist(genevec)
      genevec=gsub("mutORamp","",genevec)
      genevec=gsub("delORmut","",genevec)
      genevec=gsub("_","",genevec)
      genevec=gsub("amp","",genevec)
      genevec=gsub("mut","",genevec)
      genevec=gsub("del","",genevec)
      dfi=match(features(drugi())[[1]][[canceri]][i],results[[drugi()]][[canceri]]$df$genes)
      if(identical(input$highlight,"beta")){
        if(identical(input$druggroup,"CCLE")|input$drug %in% paperdrugs){
          value1a=results[[drugi()]][[canceri]]$df$effect[dfi]
        Fn=ecdf(abs(results[[drugi()]][[canceri]]$df$effect))
          value1=Fn(abs(value1a))
        }
        else if(identical(input$druggroup,"Sanger")){
          effecti=match(features(drugi())[[1]][[canceri]][i],names(results[[drugi()]][[canceri]]$effect))
          value1a=as.numeric(results[[drugi()]][[canceri]]$effect[effecti])
          Fn=ecdf(abs(results[[drugi()]][[canceri]]$effect))
          value1=Fn(abs(value1a))
        }
      } else if(identical(input$highlight,"freqcounts")){
        value1=results[[drugi()]][[canceri]]$df$freqCounts[dfi]
      }
      else if(identical(input$highlight,"freqevents")){
        value1=DrugCancer(canceri)$df$freqEvents[dfi]
      }
      if(isTRUE(results[[drugi()]][[canceri]]$df$posFreq[dfi]<0.5)) value1=-1*value1
      gene=c(gene,genevec)
      value=c(value,rep(value1,length(genevec)))
    }
    #pgenedata=data.frame(gene,value)
    df=matrix(nrow=0,ncol=2)
    for(i in 1:length(gene)){
      if((gene[i] %in% df[,1])==F)
        if(identical(abs(value[i]),max(abs(value[gene==gene[i]]))))
          df=rbind(df,c(gene[i],value[i]))
    }
    value=as.numeric(df[,2])
    names(value)<-df[,1]
    return(value)
  })
  
  #generate pathview object, map image file
  pathviewmake=reactive({
    if(input$makepath==0)
      return(NULL)
    isolate({
      out.suffix=paste(input$cancer[1],input$drug,sep="_")
      digits=ifelse(input$highlight=="beta",5,1)
      idinput=ifelse(input$keggmethod=="main",input$keggid,input$keggid2)  
      pv.out=pathview(gene.data=pathwaygenedata(),pathway.id=idinput,out.suffix=out.suffix,kegg.native=T,
                      gene.idtype="SYMBOL",same.layer=F,gene.annotpkg="org.Hs.eg.db",limit=list(gene=c(-1,1),cpd=1),
                      both.dirs=list(gene=T,cpd=T))
      return(pv.out)
    })
  })
  
  ##display Pathview image, stats, and data table when "Make path" button is selected
output$pathwaymap=renderUI({
  if(input$makepath==0) return(NULL)
  isolate({
    if(!is.null(pathviewmake()))
        list(imageOutput("pathway1",height="750px"),
        textOutput("fisherinfo"),
        tableOutput("fisherTable"))
         })
})

  #render pathway map image file
output$pathway1=renderImage({
  idinput=ifelse(input$keggmethod=="main",input$keggid,input$keggid2)
  filename=normalizePath(file.path(paste('/gluster/home/cnoecker/vdsapp3/hsa',idinput,'.',input$cancer[1],'_',input$drug,'.png',sep='')))
    list(src=filename,width=700)
},deleteFile=FALSE)
  
  #Fisher test info on selected pathway
output$fisherinfo=renderText({
  paste("Fisher test results for pathway '",pathwayname(),"':")
})
  
  #Fisher test statistics table
output$fisherTable=renderTable({
  fisherresult=fisherresults[[drugi()]][[match(input$cancer[1],allcancers)]][fisherresults[[drugi()]][[match(input$cancer[1],allcancers)]][,1]==input$keggid2,]
  fishertab=data.frame(Cancer=input$cancer,Drug=input$drug,KeggPath=pathwayname(),
                       OddsRatio=round(as.numeric(fisherresult[2]),digits=3),PVal=round(as.numeric(fisherresult[3]),digits=5))
},include.rownames=F,align=rep("c",6))
  
  #Pathview object data table output
  output$pathwaytable=renderTable({
    if(input$makepath==0)
      return(NULL)
    isolate({
      pathviewmake()$plot.data.gene
    })
    
  })
  
  ##output content for feature-specific info
    for(i in 1:maxfeats){
    local({ #has to do with the order things get updated
      ###Feature by cancer table
      my_i =i
      textname=paste("header",my_i,sep="")
      output[[textname]]<-renderText({
        if(paperfeatures2[genefeati()[my_i]] %in% features2(drugi())[[1]])
        paste("Feature: ",paperfeatures2[genefeati()[my_i]],sep="")
      })
      tablename=paste("table",my_i,sep="")
      output[[tablename]]<-renderTable({
        feat_t=data.frame(Cancer=allcancers,
                          coefficient=rep(NA,16),
                          posFreq=rep(NA,16),
                          freqCounts=rep(NA,16),
                          freqEvents=rep(NA,16),
                          pval=rep(NA,16),stringsAsFactors=F)
        if(identical(input$druggroup,"CCLE")){
          for(j in 1:length(allcancers)){
            if(j<=length(results[[drugi()]])){
            feat_t[j,2]=results[[drugi()]][[j]]$df$effect[match(myfeat()[my_i],results[[drugi()]][[j]]$df$genes)]
            feat_t[j,3]=results[[drugi()]][[j]]$df$posFreq[match(myfeat()[my_i],results[[drugi()]][[j]]$df$genes)]
            feat_t[j,4]=results[[drugi()]][[j]]$df$freqCounts[match(myfeat()[my_i],results[[drugi()]][[j]]$df$genes)]
            feat_t[j,5:6]=results[[drugi()]][[j]]$df[match(myfeat()[my_i],results[[drugi()]][[j]]$df$genes),7:8]
            }
          }
          feat_t=feat_t[!is.na(feat_t[,5]),]
        }
        else if(identical(input$druggroup,"Sanger")){
          for(j in 1:length(allcancers)){ 
            if(input$drug %in% paperdrugs & allcancers[j] %in% papercancers){
              feat_t[j,2]=results[[drugi()]][[j]]$df$effect[match(myfeat()[my_i],results[[drugi()]][[j]]$df$genes)]
              feat_t[j,3]=results[[drugi()]][[j]]$df$posFreq[match(myfeat()[my_i],results[[drugi()]][[j]]$df$genes)]
              feat_t[j,4]=results[[drugi()]][[j]]$df$freqCounts[match(myfeat()[my_i],results[[drugi()]][[j]]$df$genes)]
              feat_t[j,5:6]=results[[drugi()]][[j]]$df[match(myfeat()[my_i],results[[drugi()]][[j]]$df$genes),7:8]
              }else {
                if(j<=length(results[[drugi()]])){
              feat_t[j,2]=results[[drugi()]][[j]]$effect[match(myfeat()[my_i],names(results[[drugi()]][[j]]$effect))]
            feat_t[j,3]=results[[drugi()]][[j]]$df[match(myfeat()[my_i],results[[drugi()]][[j]]$df$genes),3]
            feat_t[j,4]=results[[drugi()]][[j]]$df[match(myfeat()[my_i],results[[drugi()]][[j]]$df$genes),4]
            feat_t[j,5:6]=results[[drugi()]][[j]]$df[match(myfeat()[my_i],results[[drugi()]][[j]]$df$genes),6:7]
                }
              }
          }
          feat_t=feat_t[!is.na(feat_t[,5]),]
        }
        if(dim(feat_t)[1]<1) return(NULL) else return(feat_t)
      },align=rep("c",7),digits=c(0,0,4,2,2,2,4),include.rownames=F)
      
      ###text and tables by drug for each cancer
      
      blcatextname=paste("blcafreq",my_i,sep="")
      output[[blcatextname]]<-renderText({
        if(myfeat()[my_i] %in% paperfeatures[[1]])
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[1,my_i],sep="")      
      })
      blcatablename=paste("blcatable",my_i,sep="")
      output[[blcatablename]]<-renderTable({
        if(myfeat()[my_i] %in% paperfeatures[[1]]&!is.null(input$num_drugs)){
          feattable=featbydrug(x=1,onefeat=myfeat()[my_i])
          if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      kirctextname=paste("kircfreq",my_i,sep="")
      output[[kirctextname]]<-renderText({
        if(myfeat()[my_i] %in% paperfeatures[[2]])
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[2,my_i],sep="")
      })
      kirctablename=paste("kirctable",my_i,sep="")
      output[[kirctablename]]<-renderTable({
        if(myfeat()[my_i] %in% paperfeatures[[2]]&!is.null(input$num_drugs)){
          feattable=featbydrug(x=2,onefeat=myfeat()[my_i])
        if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      gbmtextname=paste("gbmfreq",my_i,sep="")
      output[[gbmtextname]]<-renderText({
        if(myfeat()[my_i] %in% paperfeatures[[3]])
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[3,my_i],sep="")
      })
      gbmtablename=paste("gbmtable",my_i,sep="")
      output[[gbmtablename]]<-renderTable({
        if(myfeat()[my_i] %in% paperfeatures[[3]]&!is.null(input$num_drugs)){
          feattable=featbydrug(x=3,onefeat=myfeat()[my_i])
        if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      lamltextname=paste("lamlfreq",my_i,sep="")
      output[[lamltextname]]<-renderText({
        if(myfeat()[my_i] %in% paperfeatures[[4]])
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[4,my_i],sep="")
        
      })
      lamltablename=paste("lamltable",my_i,sep="")
      output[[lamltablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[4]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=4,onefeat=myfeat()[my_i])
        if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      crctextname=paste("crcfreq",my_i,sep="")
      output[[crctextname]]<-renderText({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[5]]))
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[5,my_i],sep="")
      })
      crctablename=paste("crctable",my_i,sep="")
      output[[crctablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[5]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=5,onefeat=myfeat()[my_i])
        if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
}
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      ovtextname=paste("ovfreq",my_i,sep="")
      output[[ovtextname]]<-renderText({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[6]]))
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[6,my_i],sep="")
      })
      ovtablename=paste("ovtable",my_i,sep="")
      output[[ovtablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[6]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=6,onefeat=myfeat()[my_i])
        if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,]) 
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      ucectextname=paste("ucecfreq",my_i,sep="")
      output[[ucectextname]]<-renderText({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[7]]))
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[7,my_i],sep="")
      })
      ucectablename=paste("ucectable",my_i,sep="")
      output[[ucectablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[7]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=7,onefeat=myfeat()[my_i])
        if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      pradtextname=paste("pradfreq",my_i,sep="")
      output[[pradtextname]]<-renderText({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[8]]))
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[8,my_i],sep="")
      })
      pradtablename=paste("pradtable",my_i,sep="")
      output[[pradtablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[8]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=8,onefeat=myfeat()[my_i])
          if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      skcmtextname=paste("skcmfreq",my_i,sep="")
      output[[skcmtextname]]<-renderText({
        if(myfeat()[my_i] %in% paperfeatures[[9]])
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[9,my_i],sep="")
      })
      skcmtablename=paste("skcmtable",my_i,sep="")
      output[[skcmtablename]]<-renderTable({
        if(myfeat()[my_i] %in% paperfeatures[[9]]&!is.null(input$num_drugs)){
          feattable=featbydrug(x=9,onefeat=myfeat()[my_i])
          if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
          }
        },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      lusctextname=paste("luscfreq",my_i,sep="")
      output[[lusctextname]]<-renderText({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[10]]))
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[10,my_i],sep="")
      })
      lusctablename=paste("lusctable",my_i,sep="")
      output[[lusctablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[10]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=10,onefeat=myfeat()[my_i])
          if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
          },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      luadtextname=paste("luadfreq",my_i,sep="")
      output[[luadtextname]]<-renderText({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[11]]))
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[11,my_i],sep="")
      })
      luadtablename=paste("luadtable",my_i,sep="")
      output[[luadtablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[11]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=11,onefeat=myfeat()[my_i])
          if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
          },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      brcatextname=paste("brcafreq",my_i,sep="")
      output[[brcatextname]]<-renderText({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[12]]))
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[12,my_i],sep="")
      })
      brcatablename=paste("brcatable",my_i,sep="")
      output[[brcatablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[12]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=12,onefeat=myfeat()[my_i])
          if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      brca.3negtextname=paste("brca.3negfreq",my_i,sep="")
      output[[brca.3negtextname]]<-renderText({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[13]]))
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[13,my_i],sep="")
      })
      brca.3negtablename=paste("brca.3negtable",my_i,sep="")
      output[[brca.3negtablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[13]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=13,onefeat=myfeat()[my_i])
          if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      brca.erprtextname=paste("brca.erprfreq",my_i,sep="")
      output[[brca.erprtextname]]<-renderText({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[14]]))
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[14,my_i],sep="")
      })
      brca.erprtablename=paste("brca.erprtable",my_i,sep="")
      output[[brca.erprtablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[14]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=14,onefeat=myfeat()[my_i])
          if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      brca.her2textname=paste("brca.her2freq",my_i,sep="")
      output[[brca.her2textname]]<-renderText({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[15]]))
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[15,my_i],sep="")
      })
      brca.her2tablename=paste("brca.her2table",my_i,sep="")
      output[[brca.her2tablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[15]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=15,onefeat=myfeat()[my_i])
          if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
      
      stadtextname=paste("stadfreq",my_i,sep="")
      output[[stadtextname]]<-renderText({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[16]]))
          paste(paperfeatures2[genefeati()[my_i]],"\n\n Frequency in TCGA population: ",featurefreq()[16,my_i],sep="")
      })
      stadtablename=paste("stadtable",my_i,sep="")
      output[[stadtablename]]<-renderTable({
        if(isTRUE(myfeat()[my_i] %in% paperfeatures[[16]])&!is.null(input$num_drugs)){
          feattable=featbydrug(x=16,onefeat=myfeat()[my_i])
          if(dim(feattable)[1]<input$num_drugs) return(feattable) else return(feattable[1:input$num_drugs,])
        }
      },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
    })
  
    
    #Obtain model information for one drug
    feat1drug=function(dr,canc,onefeat){
      if(length(results[[dr]])<canc) return(NULL) else {
      if(drugs2[dr] %in% ccledrugs){
        dfi=match(onefeat,results[[dr]][[canc]]$df$genes)
        return(c(Drug=drugs[dr],
                 coefficient=round(results[[dr]][[canc]]$df[dfi,4],digits=4),
                 positiveFreq=round(results[[dr]][[canc]]$df[dfi,3],digits=2),
                 freqCounts=round(results[[dr]][[canc]]$df[dfi,5],digits=4),           
                 pval=round(results[[dr]][[canc]]$df[dfi,8],digits=4)))
      } else if(drugs2[dr] %in% sangerdrugs){
        if((drugs2[dr] %in% paperdrugs)&&(allcancers[canc] %in% papercancers)){
          dfi=match(onefeat,results[[dr]][[canc]]$df$genes)
          return(c(Drug=drugs[dr],
                   coefficient=round(results[[dr]][[canc]]$df[dfi,4],digits=4),
                   positiveFreq=round(results[[dr]][[canc]]$df[dfi,3],digits=2),
                   freqCounts=round(results[[dr]][[canc]]$df[dfi,5],digits=4),           
                   pval=round(results[[dr]][[canc]]$df[dfi,8],digits=4)))
        } else{
          dfi=match(onefeat,results[[dr]][[canc]]$df$genes)
        return(c(Drug=drugs[dr],
                 coefficient=round(results[[dr]][[canc]]$effect[match(onefeat,names(results[[dr]][[canc]]$effect))],digits=4),
                 positiveFreq=round(results[[dr]][[canc]]$df[dfi,3],digits=2),
                 freqCounts=round(results[[dr]][[canc]]$df[dfi,4],digits=4),           
                 pval=round(results[[dr]][[canc]]$df[dfi,7],digits=4)))
        }
      }
      }
    }
    
    ##Obtain model information by drug
    featbydrug=function(x,onefeat){
      # feat_t_drug=data.frame(Drug=drugs,effect=rep(NA,length(drugs)),counts=rep(NA,length(drugs)),posFreq=rep(NA,length(drugs)),freqCounts=rep(NA,length(drugs)),noEvents=rep(NA,length(drugs)),freqEvents=rep(NA,length(drugs)),pval=rep(NA,length(drugs)),stringsAsFactors=F)
      if(allcancers[x] %in% papercancers[10:13]|input$paperonly){ 
        paperdrugsi=match(paperdrugs,drugs2)
        featalldrugs=sapply(paperdrugsi,feat1drug,canc=x,onefeat=onefeat)
      } else
        featalldrugs=sapply(1:length(drugs),feat1drug,canc=x,onefeat=onefeat) 
      featalldrugs=t(featalldrugs)
      if(input$drugsort=="freqcounts")
        featalldrugs=data.frame(featalldrugs[order(featalldrugs[,4],decreasing=T),]) else
          if(input$drugsort=="beta") {
            abscoef=sapply(featalldrugs[,2],as.numeric)
            abscoef=sapply(abscoef,abs)
            featalldrugs=featalldrugs[order(abscoef,decreasing=T),]} else
              if(input$drugsort=="pval")
                featalldrugs=featalldrugs[order(featalldrugs[,5],decreasing=F),]
      featalldrugs=featalldrugs[!is.na(featalldrugs[,2]),]
      return(featalldrugs)
      #  for(i in 1:length(drugs)){
      #    feat_t_drug[i,2]=results[[i]][[x]]$effect[match(input$feature,names(results[[i]][[x]]$effect))]
      #   feat_t_drug[i,3:8]=results[[i]][[x]]$df[match(input$feature,results[[i]][[x]]$df$genes),2:7]
      #  }
    }
    
    #obtain frequency of a genetic feature in TCGA population
    featurefreq=reactive({
      if(!is.null(input$gene)){
        featurefreqs=matrix(rep(NA),nrow=length(allcancers),ncol=length(genefeati()))
        for(i in 1:length(allcancers)){
          featurefreqs[i,]=sapply(1:length(genefeati()),function(x){
            round(results[[4]][[i]]$df$freqEvents[match(myfeat()[x],results[[4]][[i]]$df$genes)],digits=4)
          })}
        return(featurefreqs)
      }
    })
    
    
    
  }
  
  # cell_html <- function(table_cell) paste0('<td>', table_cell, '</td>')
  # 
  # row_html <- function(table_row) {
  #   cells <- sapply(table_row, cell_html)
  #   collapse_cells <- paste0(cells, collapse='')
  #   paste0('<tr>', collapse_cells, '</tr>')
  # }
  # 
  # output$effects1c=renderText({
  #   df_rows <- apply(effects_tab(1), 1, row_html)
  #   
  # })
  
  # output$feat_table_drug1=renderTable({
  #     if(isTRUE(input$feature %in% features[[1]]))
  #      featbydrug(1)[1:input$num_drugs,] 
  #   },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq1=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[1])
  # })
  # output$feat_table_drug2=renderTable({
  #   if(isTRUE(input$feature %in% features[[2]])){
  #     featbydrug(2)[1:input$num_drugs,] 
  #   }
  # },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq2=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[2])
  # })
  # 
  # output$feat_table_drug3=renderTable({
  #   if(isTRUE((input$feature %in% features[[3]]))){
  #     featbydrug(3)[1:input$num_drugs,] 
  #   }
  # },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq3=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[3])
  # })
  # output$feat_table_drug4=renderTable({
  #   if(isTRUE((input$feature %in% features[[4]]))){
  #     featbydrug(4)[1:input$num_drugs,] 
  #   }
  # },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq4=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[4])
  # })
  # output$feat_table_drug5=renderTable({
  #   if(isTRUE((input$feature %in% features[[5]]))){
  #     featbydrug(5)[1:input$num_drugs,] 
  #   }
  # },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq5=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[5])
  # })
  # output$feat_table_drug6=renderTable({
  #   if(isTRUE((input$feature %in% features[[6]]))){
  #     featbydrug(6)[1:input$num_drugs,] 
  #   }
  # },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq6=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[6])
  # })
  # output$feat_table_drug7=renderTable({
  #   if(isTRUE((input$feature %in% features[[7]]))){
  #     featbydrug(7)[1:input$num_drugs,] 
  #   }
  # },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq7=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[7])
  # })
  # output$feat_table_drug8=renderTable({
  #   if(isTRUE((input$feature %in% features[[8]]))){
  #     featbydrug(8)[1:input$num_drugs,] 
  #   }
  # },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq8=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[8])
  # })
  # output$feat_table_drug9=renderTable({
  #   if(isTRUE((input$feature %in% features[[9]]))){
  #     featbydrug(9)[1:input$num_drugs,] 
  #   }
  # },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq9=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[9])
  # })
  # output$feat_table_drug10=renderTable({
  #   if(isTRUE((input$feature %in% features[[10]]))){
  #     featbydrug(10)[1:input$num_drugs,] 
  #   }
  # },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq10=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[10])
  # })
  # output$feat_table_drug11=renderTable({
  #   if(isTRUE((input$feature %in% features[[11]]))){
  #     featbydrug(11)[1:input$num_drugs,] 
  #   }
  # },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq11=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[11])
  # })
  # output$feat_table_drug12=renderTable({
  #   if(isTRUE((input$feature %in% features[[12]]))){
  #     featbydrug(12)[1:input$num_drugs,] 
  #   }
  # },include.rownames=F,digits=c(0,0,4,2,2,4),align=rep("c",6))
  # output$feat_freq12=renderText({
  #   paste("Frequency of aberration in TCGA population:",featurefreq()[12])
  # })
  # 
  # output$featDrugPlots=renderPlot({
  #   if(numcanc()>0){
  #   par(mfrow=c(1,numcanc()))
  #   if(input$plotdrugs==0)
  #     return(NULL)
  #   isolate({
  #   for(i in 1:numcanc()){
  #     xvar=c()
  #     yvar=c()
  #     betasdrug=featbydrug(indices()[i])[,2]
  #     freqcountsdrug=featbydrug(indices()[i])[,4]
  #     pvalsdrug=featbydrug(indices()[i])[,5]
  #     if(input$xaxis=="beta") xvar<-betasdrug else 
  #       if(input$xaxis=="freqcounts") xvar<-freqcountsdrug else
  #         if(input$xaxis=="pval") xvar<-pvalsdrug
  #     if(input$yaxis=="beta") yvar<-betasdrug else 
  #       if(input$yaxis=="freqcounts") yvar<-freqcountsdrug else
  #         if(input$yaxis=="pval") yvar<-pvalsdrug
  #   plot(xvar[featbydrug(i)[,3]<0.3],yvar[featbydrug(i)[,3]<0.3],col=c("blue"),xlab=input$xaxis,ylab=input$yaxis,xlim=c(-1,1),ylim=c(0,1))
  #     points(xvar[featbydrug(i)[,3]>=0.3&featbydrug(i)[,3]>=0.3&featbydrug(i)[,3]<=0.7],yvar[featbydrug(i)[,3]>=0.3&featbydrug(i)[,3]<=0.7],col=c("black"))
  #     points(xvar[featbydrug(i)[,3]>0.7],yvar[featbydrug(i)[,3]>0.7],col=c("red"))
  #     legend("bottomleft",c("Negative Effect","Neutral","Positive Effect"),pch=1,col=c("blue","black","red"))
  #     }
  # })
  #   }
  #   })
  
})
#output$auc_celllines=renderText(vdsObj$auc_celllines)
#output$metrics=renderText(vdsObj$diseaseFMAT$crc$metric)
#output$perf_metrics=renderTable(
#data.frame(Value=c(vdsObj$auc_celllines,vdsObj$diseaseFMAT$crc$metric[1],vdsObj$diseaseFMAT$crc$metric[2],vdsObj$diseaseFMAT$crc$N),row.names=c("AUC all cell lines","AUC","Spearman coefficient","N TCGA samples")))

#})

