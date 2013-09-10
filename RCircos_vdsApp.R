##### RCircos code for Shiny app #####
library(RCircos)

#RCircos draw chromosome function-adjusted to speed up, remove extra stuff

RCircos.Chromosome.Ideogram.Plot=function () 
{
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  right.side <- nrow(RCircos.Pos)/2
  outer.location <- RCircos.Par$chr.ideog.pos + RCircos.Par$chrom.width
  inner.location <- RCircos.Par$chr.ideog.pos
  chroms <- unique(RCircos.Cyto$Chromosome)
  for (a.chr in 1:length(chroms)) {
    the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr], 
                            ]
    start <- the.chr$Location[1] - the.chr$Unit[1] + 1
    end <- the.chr$Location[nrow(the.chr)]
    mid <- round((end - start + 1)/2, digits = 0) + start
    chr.color <- the.chr$ChrColor[nrow(the.chr)]
    pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, 
               RCircos.Pos[end:start, 1] * inner.location)
    pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, 
               RCircos.Pos[end:start, 2] * inner.location)
    polygon(pos.x, pos.y)
    chr.name <- sub(pattern = "chr", replacement = "", chroms[a.chr])
    text(RCircos.Pos[mid, 1] * RCircos.Par$chr.name.pos, 
         RCircos.Pos[mid, 2] * RCircos.Par$chr.name.pos, label = chr.name, 
         srt = RCircos.Pos$degree[mid])
    lines(RCircos.Pos[start:end, ] * RCircos.Par$highlight.pos, 
          col = chr.color, lwd = RCircos.Par$highlight.width)
  }
for (a.band in 1:nrow(RCircos.Cyto)) {
a.color <- RCircos.Cyto$BandColor[a.band]
 if (a.color == "white") {
   next
 }
 start <- RCircos.Cyto$Location[a.band] - RCircos.Cyto$Unit[a.band] + 
   1
 end <- RCircos.Cyto$Location[a.band]
 pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, 
            RCircos.Pos[end:start, 1] * inner.location)
 pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, 
            RCircos.Pos[end:start, 2] * inner.location)
 polygon(pos.x, pos.y, col = a.color, border = NA)
}
}

data(UCSC.HG19.Human.CytoBandIdeogram)
data(RCircos.Tile.Data)
data(RCircos.Gene.Label.Data)
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
gene.data <- RCircos.Gene.Label.Data
RCircos.Set.Core.Components(cyto.info,chr.exclude=NULL,tracks.inside=5,tracks.outside=0)


########  Put together complete location data for genes included in VDS models
library(org.Hs.eg.db)
complete.gene.data=select(org.Hs.eg.db,keys=genes,cols=c("CHR","CHRLOC","CHRLOCEND"),keytype="SYMBOL")
complete.gene.data=data.frame(complete.gene.data)
complete.gene.data=subset(complete.gene.data, !duplicated(SYMBOL))
complete.gene.data$CHRLOC=lapply(complete.gene.data$CHRLOC,function(x){x=abs(x)})
complete.gene.data$CHRLOCEND=lapply(complete.gene.data$CHRLOCEND,function(x){x=abs(x)})
complete.gene.data=data.frame(cbind(Chromosome=complete.gene.data$CHR,chromStart=complete.gene.data$CHRLOC,chromEnd=complete.gene.data$CHRLOCEND,Gene=complete.gene.data$SYMBOL),stringsAsFactors=F)
chromStart=as.vector(unlist(complete.gene.data$CHRLOC))
Chromosome=as.vector(unlist(complete.gene.data$CHR))
chromEnd=as.vector(unlist(complete.gene.data$CHRLOCEND))
Gene=as.vector(unlist(complete.gene.data$SYMBOL))
complete.gene.data=data.frame(Chromosome,chromStart,chromEnd,Gene,stringsAsFactors=F)
write.csv(complete.gene.data,"complete_gene_data.csv",row.names=F)

complete.gene.data=read.csv("/gluster/home/cnoecker/complete_gene_data.csv")


################## ui.R RCircos code

##Inputs for Circos Plot tab
conditionalPanel(condition="input.tabset1=='Circos Plot'",
                selectInput("filtertype","Include genes based on counts or p-value?",choices=c("Counts/Importance score","P-value"))
                ),
conditionalPanel(condition="input.tabset1=='Circos Plot' && input.filtertype=='P-Value'",
                numericInput("pthreshold","Include features with what level of significance?",value=0.05,min=0,max=1,step=0.01)),
conditionalPanel(condition="input.tabset1=='Circos Plot' && input.filtertype=='Counts/Importance score'",
               numericInput("cthreshold","Minimum importance score of features included?",value=0.3,min=0,max=1,step=0.05)),
conditionalPanel(condition="input.tabset1=='Circos Plot'",
              actionButton("makecircos","Make Circos Plot(s)")),

#display plot
tabPanel("Circos Plot",
        plotOutput("CircosPlot",height="600px")),

################### server.R RCircos code

##Circos Plot Outputs
circplot=function(x){
  new.gene.data1<-data.frame(Chromosome=rep("",0),chromStart=rep("",0),chromEnd=rep("",0),Gene=rep("",0),stringsAsFactors=F)
    for(i in 1:length(complete.gene.data$Gene)){
    g=grep(complete.gene.data$Gene[i],DrugCancer(indices()[x])$df$genes)
    if(length(g)>0){
      for(j in 1:length(g)){
        if(identical(input$filtertype,"Counts/Importance score")){
          if(DrugCancer(indices()[x])$df$freqCounts[g[j]]>input$cthreshold)
            new.gene.data1=rbind(new.gene.data1,complete.gene.data[i,]) } else 
        if(identical(input$filtertype,"P-Value")){
            if(DrugCancer(indices()[x])$df$pvals[g[j]]<input$pthreshold)
              new.gene.data1=rbind(new.gene.data1,complete.gene.data[i,])
            }
      }
    }
  }
  new.gene.data1=subset(new.gene.data1,!duplicated(new.gene.data1$Gene))
  #   for(i in 1:length(complete.gene.data$Gene)){
  #     if(!is.na(pmatch(complete.gene.data$Gene[i],features[[match(input$cancer,cancers)]],nomatch=NA)))
  #       if(DrugCancer()$df$pvals[match(complete.gene.data$Gene[i],DrugCancer()$df$genes)]<0.05)
  #       new.gene.data=rbind(new.gene.data,complete.gene.data[i,])
  #   }
  
}

output$CircosPlot=renderPlot({
    if(input$makecircos==0)
      return(NULL)
    isolate({
      new.gene.data=do.call(rbind,lapply(1:numcanc(),circplot))
      RCircos.Set.Plot.Area()
      RCircos.Chromosome.Ideogram.Plot()
      RCircos.Gene.Name.Plot(new.gene.data,name.col=4,track.num=2,side="in")
        })
 
},height="auto")


