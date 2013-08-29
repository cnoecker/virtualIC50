
mafDir="/home/jguinney/data/tcga_mafs/"
mafFiles <- list(
                 #blca="BLCA-28-original.aggregated.tcga.somatic.maf",
                 #kirc="BI_and_BCM_1.4.aggregated.tcga.somatic.maf",
                 #gbm="gbm_liftover.aggregated.capture.tcga.uuid.somatic.maf", 
                 #laml="genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2.2.12.0.somatic.maf",
                 #ucec="UCEC_somatic.maf",
                 #prad="PR-TCGA-Analysis_set.aggregated.tcga.somatic.maf",
                 skcm="TCGA_SKCM_266_PAIR.aggregated.capture.tcga.uuid.somatic.maf",
                 #lusc="LUSC_Paper_v8.aggregated.tcga.somatic.maf",
                 luad="PR_TCGA_LUAD_PAIR_Capture_All_Pairs_QCPASS.aggregated.capture.tcga.uuid.somatic.maf",
                 stad="PR_TCGA_STAD_PAIR_Capture_All_Pairs_QCPASS_v4.aggregated.capture.tcga.uuid.somatic.maf"
                 #brca="genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.5.3.0.somatic.maf")
)
buildMutationMatrixFromPANCAN <- function(tumorType, cancer.genes=TRUE, missenseFilter=TRUE){
  
  if(class(missenseFilter)=="list"){
    cosmicProteinPositions = missenseFilter
    missenseFilter = TRUE
  }else if(class(missenseFilter) == "logical" & missenseFilter){
    cat("loading cosmic database for missense filtering...")
    cosmicProteinPositions <- processCosmicMutationFile()
    cat("done\n")
  }
  if(tumorType %in% names(mafFiles)){
    mafFile <- paste(mafDir,mafFiles[[tumorType]],sep="")
  }else{
    e <- loadEntity("syn1710680")
    mafFile <- paste(e$cacheDir,"/somatic_mafs_cleaned/",tumorType,"_cleaned.maf",sep="")
  }
  tbl <- read.table(mafFile,header=T,as.is=T,quote="",comment="#",sep="\t")
  
  pat.ids <- patient.ids(gsub("-",".",tbl$Tumor_Sample_Barcode))
  tbl$pat.ids <- pat.ids
  
  if(cancer.genes){
    cbio <- read.table("./resources/cbio_cancer_genes.txt",sep="\t",header=T,as.is=T,quote="")$Gene.Symbol
    sanger <- read.table("./resources/cancer_gene_census.txt",sep="\t",header=T,as.is=T,quote="")$Symbol
    genes <- union(union(cbio, sanger),"MTOR")
  }else{
    genes <- sort(unique(tbl$Hugo_Symbol))
  }
  
  filter.mask <- tbl$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation") & tbl$Hugo_Symbol %in% genes
  tblheader <- colnames(tbl)
  aachange <- ifelse("Protein_Change" %in% tblheader,"Protein_Change",
                     tblheader[grep("amino_acid_change", tblheader)])
  filtered.tbl <- tbl[filter.mask,c("Tumor_Sample_Barcode","pat.ids","Hugo_Symbol","Variant_Classification",aachange)]
  mafGenes <- sort(unique(filtered.tbl$Hugo_Symbol))
  genes <- intersect(mafGenes, genes)
  samples <- sort(unique(filtered.tbl$pat.ids))
  
  M <- matrix("", nrow=length(genes),ncol=length(samples),dimnames=list(genes,samples))
  
  idxs.list <- oneToManyOperation(genes, filtered.tbl$Hugo_Symbol)
  
  for(i in 1:length(idxs.list)){
    gene <- genes[[i]]
    filteredMAFByGene <- filtered.tbl[idxs.list[[i]],]
    
    if(missenseFilter & length(cosmicProteinPositions[[gene]]) > 0){
      posCounts <- table(cosmicProteinPositions[[gene]])
      x <- rep(0, max(as.numeric(names(posCounts))))
      x[as.numeric(names(posCounts))] <- posCounts
      # sliding count window of size 10, so both 4 AA pos on both side
      y <- runsum(x, 9)
      
      missense.mask <- filteredMAFByGene$Variant_Classification == "Missense_Mutation"
      aa_pos <- as.numeric(gsub("p\\.\\w(\\d+).*","\\1",filteredMAFByGene[missense.mask,aachange])) - 4
      aa_pos[aa_pos < 1] <- 1
      tmpmask <- rep(TRUE, sum(missense.mask))
      tmpmask[aa_pos > length(y)] <- FALSE
      # filter out if window has 2 or fewer counts
      tmpmask[y[aa_pos] <= 2] <- FALSE
      keepmissense.mask <- missense.mask
      keepmissense.mask[keepmissense.mask] <- tmpmask
      
      filteredMAFByGene <- filteredMAFByGene[!missense.mask | keepmissense.mask, ]
    }
    
    if(nrow(filteredMAFByGene) > 0){
      M[i,] <- unlist(oneToManyOperation(samples,filteredMAFByGene$pat.ids, function(idxs){
        ifelse(is.null(idxs),"",paste(sort(unique(filteredMAFByGene[idxs,aachange])),collapse=";"))
      }))
    }
  }
  return (M)
}
  
###############################################
###  
  
buildMutationMatrix <- function(synapseId, cancer.genes=TRUE){
  if(startsWith(synapseId,"syn")){  
    e <- loadEntity(synapseId)
    synapseId <- paste(e$cacheDir,"/",e$files,sep="")
  }
  tbl <- read.table(synapseId,header=T,as.is=T,quote="",comment="",sep="\t")
  
  pat.ids <- patient.ids(gsub("-",".",tbl$Tumor_Sample_Barcode))
  tbl$pat.ids <- pat.ids
  
  if(cancer.genes){
    cbio <- read.table("./resources/cbio_cancer_genes.txt",sep="\t",header=T,as.is=T,quote="")$Gene.Symbol
    sanger <- read.table("./resources/cancer_gene_census.txt",sep="\t",header=T,as.is=T,quote="")$Symbol
    genes <- union(union(cbio, sanger),"MTOR")
  }else{
    genes <- sort(unique(tbl$Hugo_Symbol))
  }
  
  filter.mask <- tbl$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation") & tbl$Hugo_Symbol %in% genes
  protein_change_header <- ifelse("amino_acid_change_WU" %in% colnames(tbl), "amino_acid_change_WU", 
                                  ifelse("AAChange" %in% colnames(tbl), "AAChange", "Protein_Change"))
  filtered.tbl <- tbl[filter.mask,c("Tumor_Sample_Barcode","pat.ids","Hugo_Symbol","Variant_Classification",protein_change_header)]
  mafGenes <- sort(unique(filtered.tbl$Hugo_Symbol))
  genes <- intersect(mafGenes, genes)
  samples <- sort(unique(filtered.tbl$pat.ids))
  
  
  M <- matrix("", nrow=length(genes),ncol=length(samples),dimnames=list(genes,samples))
  
  idxs.list <- oneToManyOperation(genes, filtered.tbl$Hugo_Symbol)
  
  for(i in 1:length(idxs.list)){
    gene <- genes[[i]]
    filteredMAFByGene <- filtered.tbl[idxs.list[[i]],]
    
    M[i,] <- unlist(oneToManyOperation(samples,filteredMAFByGene$pat.ids, function(idxs){
      ifelse(is.null(idxs),"",paste(sort(unique(filteredMAFByGene[idxs,protein_change_header])),collapse=";"))
    }))
  }
  return (M)
}

## load cosmic mutation database and provide counts (per gene) for each position
processCosmicMutationFile <- function(file="./data~/cosmic/CosmicMutantExport_v65_280513.tsv"){
  tbl <- read.table("data~/cosmic/CosmicMutantExport_v65_280513.tsv",header=T,sep="\t",as.is=T,quote="")
  mtbl <- tbl[grepl("Missense", tbl$Mutation.Description),]
  genes <- unique(mtbl$Gene.name)
  proteinPositions <- oneToManyOperation(genes, mtbl$Gene.name, function(idxs){
    y <- mtbl[idxs,"Mutation.AA"]
    pos <- sort(as.numeric(gsub("p\\.\\w(\\d+).*","\\1",y)))
    pos
  })
  names(proteinPositions) <- genes
  return (proteinPositions)
}

# running window of counts for window size k
runsum = function(x, k) { 
  n = length(x) 
  y = x[ k:n ] - x[ c(1,1:(n-k)) ] # this is a difference from the previous cell 
  y[1] = sum(x[1:k]); # find the first sum 
  y = cumsum(y) # apply precomputed differences 
  return(y) 
}

process_firehose_mutation_package <- function(){
  mafdir <- "./data~/crc/gdac.broadinstitute.org_COADREAD.Mutation_Packager_Calls.Level_3.2013052300.0.0"
  setwd(mafdir)
  
  unlink("CRC_maf.txt")
  mafs <- dir(".",pattern="TCGA.*maf\\.txt")
  system(paste("cp ", mafs[1], " CRC_maf.txt"))
  for(i in 2:length(mafs)){
    system(paste("sed -e '1d' ", mafs[i], " >> CRC_maf.txt",sep=""))
  }
}  
