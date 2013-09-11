
############### Experimental/incomplete functions

#### server.R
### make rows of table different colors
cell_html <- function(table_cell) paste0('<td>', table_cell, '</td>')

row_html <- function(table_row) {
  cells <- sapply(table_row, cell_html)
  collapse_cells <- paste0(cells, collapse='')
  paste0('<tr>', collapse_cells, '</tr>')
}

output$effects1c=renderText({
  df_rows <- apply(effects_tab(1), 1, row_html)
  
})

######plot results for different drugs on same plot
output$featDrugPlots=renderPlot({
  if(numcanc()>0){
  par(mfrow=c(1,numcanc()))
  if(input$plotdrugs==0)
    return(NULL)
  isolate({
  for(i in 1:numcanc()){
    xvar=c()
    yvar=c()
    betasdrug=featbydrug(indices()[i])[,2]
    freqcountsdrug=featbydrug(indices()[i])[,4]
    pvalsdrug=featbydrug(indices()[i])[,5]
    if(input$xaxis=="beta") xvar<-betasdrug else 
      if(input$xaxis=="freqcounts") xvar<-freqcountsdrug else
        if(input$xaxis=="pval") xvar<-pvalsdrug
    if(input$yaxis=="beta") yvar<-betasdrug else 
      if(input$yaxis=="freqcounts") yvar<-freqcountsdrug else
        if(input$yaxis=="pval") yvar<-pvalsdrug
  plot(xvar[featbydrug(i)[,3]<0.3],yvar[featbydrug(i)[,3]<0.3],col=c("blue"),xlab=input$xaxis,ylab=input$yaxis,xlim=c(-1,1),ylim=c(0,1))
    points(xvar[featbydrug(i)[,3]>=0.3&featbydrug(i)[,3]>=0.3&featbydrug(i)[,3]<=0.7],yvar[featbydrug(i)[,3]>=0.3&featbydrug(i)[,3]<=0.7],col=c("black"))
    points(xvar[featbydrug(i)[,3]>0.7],yvar[featbydrug(i)[,3]>0.7],col=c("red"))
    legend("bottomleft",c("Negative Effect","Neutral","Positive Effect"),pch=1,col=c("blue","black","red"))
    }
})
  }
  })


############## Show results for drugs with the same target - never fully implemented

targetinfo=read.csv("/gluster/home/cnoecker/drugtargets.csv")
head(targetinfo)
table(targetinfo$Target)
drugs2=gsub("ccle_","",drugs)
drugs2=gsub("sanger_","",drugs2)

targets=rep(NA,length(drugs2))
for(i in 1:length(drugs2)){
  ind=match(drugs2[i],targetinfo$Drug)
    targets[i]=as.character(targetinfo$Target[ind])
}

write.csv(cbind(drugs2,targets),file="/gluster/home/cnoecker/targets.csv")


targetinfo=read.csv("/gluster/home/cnoecker/targets.csv")
targets=targetinfo$targets
targets=gsub(",","",targets)


##################### ui.R
#show results for drugs with same target
conditionalPanel(condition="input.tabset1=='Compare Feature Across Cancers'",
                actionButton("reldrugs","See results for this gene for other related drugs")),

#plotting options for results for different drugs
conditionalPanel(condition="input.tabset1=='Compare Feature Across Drugs'",
               selectInput("xaxis","Plot X variable",choices=c("Effect size"="beta","Importance score"="freqcounts","P-value"="pval")),
              selectInput("yaxis","Plot Y variable",choices=c("Effect size"="beta","Importance score"="freqcounts","P-value"="pval")),
             actionButton("plotdrugs","Plot Drugs")),
