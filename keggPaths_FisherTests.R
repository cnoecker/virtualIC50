############## Calculate the most enriched KEGG pathways for each cancer/drug combination

####### Get info from KEGGREST package on genes in each pathway

##access info on kegg pathways
KEGG map genes - kegg IDs, genes in each map
keggpaths=keggList("pathway","hsa") #275
keggids=names(keggpaths)
keggids=gsub("path:hsa","",keggids) #Just ID numbers

names(keggpaths)=keggids
keggpaths2=keggids
names(keggpaths2)=keggpaths
allkegggenes=keggList("hsa") #26269

kegggenes=keggLink("pathway","hsa") #22143

kegggenenames=names(kegggenes)
kegggenenames=gsub("hsa:","",kegggenenames)
genesbypath=table(kegggenes)
names(genesbypath)=gsub("path:hsa","",names(genesbypath))

x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])

kegg2=c() ##unique gene names in pathways
for(j in 1:length(kegggenenames)){
kegg2[j]=xx[names(xx)==kegggenenames[j]][[1]]
}
kegg2=unique(kegg2) ##6716

pathwaygenes=cbind(kegggenenames,kegggenes)
pathwaygenes[,2]=gsub("path:hsa","",pathwaygenes[,2])
for(j in 1:length(pathwaygenes[,1])){
  pathwaygenes[j,1]=xx[names(xx)==pathwaygenes[j,1]][[1]]
}
save(keggpaths,keggpaths2,keggids,genesbypath,kegg2,pathwaygenes,fisherresults,file="/gluster/home/cnoecker/keggpathways.rda")


####### Function to perform Fisher test for one model and one pathway

fishertest=function(keggi,drugi,canceri){
  modgenes=c()
  if(!is.na(match(drugs2[drugi],paperdrugs))&
       !is.na(match(allcancers[canceri],papercancers))){
    modgenes=unlist(sapply(1:length(paperfeatures[[canceri]]),function(i){
      dfi=match(paperfeatures[[canceri]][i],results[[drugi]][[canceri]]$df$genes)
      if(abs(results[[drugi]][[canceri]]$df$effect[dfi])>0)
        return(row.names(results[[drugi]][[canceri]]$df[dfi,])) else return(NULL)
    }))
  } else 
    if(drugi<25){
      modgenes=unlist(sapply(1:length(cclefeatures[[canceri]]),function(i){
        dfi=match(cclefeatures[[canceri]][i],results[[drugi]][[canceri]]$df$genes)
        if(abs(results[[drugi]][[canceri]]$df$effect[dfi])>0)
          return(row.names(results[[drugi]][[canceri]]$df[dfi,])) else return(NULL)
      }))
    } else if(drugi>=25){
      modgenes=unlist(sapply(1:length(sangerfeatures[[canceri]]),function(i){
        dfi=match(sangerfeatures[[canceri]][i],results[[drugi]][[canceri]]$df$genes)
        effect1=results[[drugi]][[canceri]]$effect[match(results[[drugi]][[canceri]]$df$genes[dfi],names(results[[drugi]][[canceri]]$effect))]
        if(effect1!=0)
          return(row.names(results[[drugi]][[canceri]]$df[dfi,])) else return(NULL)
      }))
    }
  pathgenes=pathwaygenes[,1][pathwaygenes[,2]==keggids[keggi]]
  modpathgenes=c()
  modpathgenes=unlist(sapply(1:length(pathgenes),function(j){
    genei=c(grep(paste("^",pathgenes[j],"_",sep=""),modgenes),grep(paste(":",pathgenes[j],"_",sep=""),modgenes))
    if(length(genei)>0) return(pathgenes[j]) else return(NULL)
  }))
  fisher.table=cbind(c(length(modpathgenes),length(modgenes)-length(modpathgenes)),c(length(pathgenes),length(kegg2)-length(pathgenes)))
  fisher1=fisher.test(fisher.table)
  return(c(fisher1$estimate,fisher1$p.value))
}



########Loop over results for all drugs

for(i in 1:length(results)){
  drugresults=list(NULL)  ### list of results for 1 drug, all cancers
  for(j in 1:length(results[[i]])){
    oddsratios=sapply(1:length(keggids),fishertest,drugi=i,canceri=j)
    oddsratios=t(oddsratios)
    ortable=cbind(keggID=keggids,oddsRatio=oddsratios[,1],pValue=oddsratios[,2])
    ortable=ortable[order(as.numeric(ortable[,2]),decreasing=T),]
    drugresults[[j]]=ortable
  }
  fisherresults[[i]]=drugresults
}

save(fisherresults,file="/gluster/home/cnoecker/fisherpathresults.rda")



