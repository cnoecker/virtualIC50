########### Read in model results

results=list(NULL)

######### Sanger results from all_2013_06_21 and CCLE from all_2013_08_03
for(i in 1:length(ccledrugs)){
  load(paste("/gluster/home/jguinney/projects/virtualIC50/results/all_2013_08_03/ccle_",ccledrugs[i],"/output.rda",sep=""))
  results[[i]]=vdsObj$diseaseFMAT
}

for(i in 1:length(sangerdrugs)){
  load(paste("/gluster/home/jguinney/projects/virtualIC50/results/all_2013_06_21/sanger_",sangerdrugs[i],"/output.rda",sep=""))
  auc_allcelllines[i+24]=vdsObj$auc_celllines
  results[[i+24]]=vdsObj$diseaseFMAT
}

############## Replace a subset with results from paper_2013_08_30
for(k in 1:length(paperdrugs)){
  drugi=grep(paperdrugs[k],drugs)[1]
  for(j in 1:length(papercancers)){
  load(paste("/gluster/home/jguinney/projects/virtualIC50/results/paper_2013_08_30/rda/",paperdrugs[k],"_",papercancers[j],".rda",sep=""))
  canceri=match(papercancers[j],cancers)
  if(is.na(canceri)) canceri=j+3
  results[[drugi]][[canceri]]=fmat
  names(results[[drugi]])[canceri]=papercancers[j]
}
}

############## Additional model performance data for paper drugs
paperperftable=data.frame(drug=paperdrugs,auc=rep(NA,length(paperdrugs)),rho=rep(NA,length(paperdrugs)),npv=rep(NA,length(paperdrugs)),
                           ppv=rep(NA,length(paperdrugs)))
for(j in 1:length(paperdrugs)){
load(paste("/gluster/home/jguinney/projects/virtualIC50/results/paper_2013_08_30/rda/",paperdrugs[j],".rda",sep=""))
  paperperftable[j,2]=vds$perf$ALL$auc
  paperperftable[j,3]=vds$perf$ALL$rho$estimate
  paperperftable[j,4]=vds$perf$ALL$npv
  paperperftable[j,5]=vds$perf$ALL$ppv
}

save(results,paperperftable,file="/gluster/home/cnoecker/vdsResults.rda")
