#Takes a Seurat Object and performs differential expression analysis 
#For each cluster/cell type
#By splitting the object into seperate samples (i.e organoids) and using DESeq2,
#which accounts for noise between samples and then looks for DEGs that overcome that noise

library(Seurat)
library(DESeq2)
library(writexl)

#Load this function as-is for use later
combineDE<-function(seur,id,condition="treat",base="wt",combineOn="orig.ident", #These are defaults, but will be overrided by what you put below
                    #Below are default settings you can change if you know what you're doing
                    minCells=20, #Minimum # of cells that must be in this cluster per sample to keep that sample
                    minBatches=2, #Minimum # of samples you can have per condition. Cannot be lower than 2.
                    minReads=10, #Mininum # of reads to have total per gene to calculate DE for that gene
                    genes=c(),  #Genes to consider for DE analysis, if you don't want to use all expressed genes.
                    form="" #design formula for DESeq2 if you want to include more variables in addition to "condition"
)
{
  print("Subsample")
  seur<-subset(seur,idents=c(id))
  
  Idents(seur)=combineOn
  genes.use=rownames(GetAssayData(seur,slot="counts"))
  if(length(genes)>0){genes.use=genes}
  
  print("Combine data per sample")
  data.all=data.frame(row.names = genes.use)
  for(i in levels(Idents(seur))) {
    temp.cells=WhichCells(seur,ident=i)
    if (length(temp.cells)==1) data.temp=(GetAssayData(seur,slot="counts")[genes.use,temp.cells])
    if (length(temp.cells)>1) data.temp=apply(GetAssayData(seur,slot="counts")[genes.use,temp.cells],1,sum)
    data.all=cbind(data.all,data.temp)
    colnames(data.all)[ncol(data.all)]=i
  }
  
  print("Filter samples for minimum cells")
  keepOrgs=names(summary(Idents(seur)))[summary(Idents(seur))>minCells]
  numOrg=length(keepOrgs)
  print(paste("Keeping", numOrg, "samples"))
  data.all=data.all[,keepOrgs]
  
  extraColumns<-strsplit(form,"+",fixed=T)[[1]]
  val=seur@meta.data[,c(condition,combineOn,extraColumns)]
  val=val[!duplicated(val[,2]),]
  rownames(val)=val[,2]
  keepBatch=as.character(val[keepOrgs,1])
  levels = levels(factor(keepBatch))
  if(length(levels)<2) {
    print("Not enough batches per treatment group with minimum # of cells!")
    return(NULL)
  }
  for (level in levels) {
    if(sum(keepBatch==level)<2) {
      print("Not enough batches per treatment group with minimum # of cells!")
      return(NULL)
    }
  }
  
  print("Save meta data")
  colDat=factor(keepBatch)
  if (base != "") { colDat = relevel(colDat, ref=base)}
  colDat=data.frame(colDat)
  colnames(colDat)="condition"
  rownames(colDat)=colnames(data.all)
  colDat[keepOrgs,extraColumns]=val[keepOrgs,extraColumns]
  
  print("Run DESeq2")
  design= ~ condition
  if(nchar(form)>0){
    design=as.formula(paste("~",form,"+ condition",sep=""))
  }
  print(design)
  
  dds <- DESeqDataSetFromMatrix(countData = data.all,colData = colDat,design = design)
  dds <- dds[ rowSums(counts(dds) > minReads)>=2, ]
  dds <- DESeq(dds)
  out=data.frame(results(dds))
  out=out[order(out$pvalue),]
  return(out)
}

#Load Seurat Object
seur = readRDS("~/Documents/DorsalKadoshima/PTEN_mito210_1m/noWT2/clusteredSeur.rds")

#Set Ident to cluster or celltypes, whatever groups you want to seperate before detecting DEGs in that group
Idents(seur) = "mergedClusts_PC30"

#Set condition to the metadata column you want DEGs between, and base to the "wildtype" or base level of that column
condition = "treat"
base = "wt"

#set combineOn to the metadata column that contains the samples (i.e. different organoids)
combineOn = "org"

#This loop will run DE analysis for each cluster and save a .xlsx file for each!
for (id in levels(seur@active.ident)) {
  print(id)
  degs <- combineDE(seur, id=id, condition=condition, base=base, combineOn=combineOn)
  degs$gene = rownames(degs)
  if (length(degs)>0) {
    write_xlsx(degs, path=paste0("clus",id,".DEGs.xlsx"), format_headers = T)
  }
}

