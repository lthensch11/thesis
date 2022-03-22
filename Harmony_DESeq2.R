library(Seurat)
library(DESeq2)
library(writexl)

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


setwd("~/Documents/Harvard/Arlotta Lab/Mito210_PTEN_1M_Harmony/")
#DEGs by cell type
Idents(PTEN.1M.combined.obj) = "CellType"
id = "CellType"
condition = "treat"
base = "wt"
combineOn = "orig.ident"

for (id in levels(PTEN.1M.combined.obj@active.ident)) {
  print(id)
  degs <- combineDE(PTEN.1M.combined.obj, id=id, condition=condition, base=base, combineOn=combineOn)
  degs$gene = rownames(degs)
  if (length(degs)>0) {
    write_xlsx(degs, path=paste0(id,".DEGs.xlsx"), format_headers = T)
  }
}

setwd("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_1M_Harmony/")
#DEGs by cell type
Idents(ARID1b.1M.combined.obj) = "CellType"
id = "CellType"
condition = "treat"
base = "wt"
combineOn = "orig.ident"

for (id in levels(ARID1b.1M.combined.obj@active.ident)) {
  print(id)
  degs <- combineDE(ARID1b.1M.combined.obj, id=id, condition=condition, base=base, combineOn=combineOn)
  degs$gene = rownames(degs)
  if (length(degs)>0) {
    write_xlsx(degs, path=paste0(id,".DEGs.xlsx"), format_headers = T)
  }
}
