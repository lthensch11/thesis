##COMPARE DEGS ACROSS CELL CLUSTERS OF INTEREST BETWEEN ASD GENE ORGS## in this case it is the pseudoSUB group of IPC, Newborn PN and Newborn DL PN

library(Seurat)
library(DESeq2)
library(writexl)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
library(viridis)


#Function that performs the DESeq calculation by summing counts in each sample
combineDE<-function(seur,id,condition="treatment",base="wt",combineOn="orig.ident",
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


setwd("~/Documents/Harvard/Arlotta Lab")

#PTEN
setwd("~/Documents/Harvard/Arlotta Lab/Mito210_PTEN_1M_rep1")
PTEN.rep1.obj <- readRDS(file = "PTEN_rep1_obj.rds")

PTEN.rep1.obj$pseudoSub = NA
PTEN.rep1.obj$pseudoSub[PTEN.rep1.obj$CellType %in% c("aRG", "Cycling Progenitors", "Subcortical", 
                                                      "Coritcal Hem", "Unknown", "Cajal-Retzius")] = FALSE
PTEN.rep1.obj$pseudoSub[PTEN.rep1.obj$CellType %in% c("Newborn PNs", "IPC", "Newborn DL PNs")] = TRUE

saveRDS(PTEN.rep1.obj, file = "celltypesSeur.rds")

setwd("~/Documents/Harvard/Arlotta Lab/Mito210_PTEN_1M_rep2")
PTEN.rep2.obj <- readRDS(file = "PTEN_rep2_obj.rds")

PTEN.rep2.obj$pseudoSub = NA
PTEN.rep2.obj$pseudoSub[PTEN.rep2.obj$CellType %in% c("aRG", "Cycling Progenitors", "Subcortical", 
                                                      "Coritcal Hem", "Unknown", "GABAergic Nuerons 1")] = FALSE
PTEN.rep2.obj$pseudoSub[PTEN.rep2.obj$CellType %in% c("Newborn PNs", "IPC", "Newborn DL PNs")] = TRUE

saveRDS(PTEN.rep2.obj, file = "celltypesSeur.rds")
saveRDS(PTEN.rep2.obj, file = "PTEN_rep2_obj.rds")


setwd("~/Documents/Harvard/Arlotta Lab/Mito210_PTEN_1M_Harmony")
PTEN_1M_harmony <- readRDS(file = "celltypesSeur.rds")

PTEN_1M_harmony$treat = NA
#create a new column called "treatment" or whatever you want and set everything to NA at first
PTEN_1M_harmony$treat[PTEN_1M_harmony$orig.ident %in% c(1,2,3)] = "mt"
PTEN_1M_harmony$treat[PTEN_1M_harmony$orig.ident %in% c(4,5)] = "wt" 

PTEN_1M_harmony$pseudoSub = NA
PTEN_1M_harmony$pseudoSub[PTEN_1M_harmony$CellType %in% c("aRG", "Cycling Progenitors", "Subcortical", 
                                                          "Coritcal Hem", "Unknown", "GABAergic Nuerons 1", "Cajal-Retzius")] = FALSE
PTEN_1M_harmony$pseudoSub[PTEN_1M_harmony$CellType %in% c("Newborn PNs", "IPC", "Newborn DL PNs")] = TRUE

saveRDS(PTEN_1M_harmony, file = "celltypesSeur.rds")
saveRDS(PTEN_1M_harmony, file = "Mito210_PTEN_1M_harmonizedObj.rds")


#ARID1b
setwd("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_1M_scRNAseq_rep1_3.15.21")
ARID1b_1M_rep1 <- readRDS(file = "clusteredSeur_ARID1B_1M.rds")

ARID1b_1M_rep1$pseudoSub = NA
ARID1b_1M_rep1$pseudoSub[ARID1b_1M_rep1$CellType %in% c("aRG", "Cycling Progenitors", "Subcortical",
                                                        "Cortical Hem", "Cajal-Retzius", "Cycling GABAergic Progenitors 1",
                                                        "GABAergic Neurons 1", "GABAergic Progenitors 1")] = FALSE
ARID1b_1M_rep1$pseudoSub[ARID1b_1M_rep1$CellType %in% c("Newborn PNs", "IPC", "Newborn DL PNs")] = TRUE

saveRDS(ARID1b_1M_rep1, file = "celltypesSeur.rds")


setwd("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_rep2_6.9.21")
Mito210_ARID1b_rep2 <- readRDS(file = "clusteredSeur.rds")

Mito210_ARID1b_rep2$pseudoSub = NA
Mito210_ARID1b_rep2$pseudoSub[Mito210_ARID1b_rep2$CellType %in% c("aRG", "Cycling Progenitors", "Subcortical",
                                                        "Cortical Hem", "Cycling GABAergic Progenitors 1",
                                                        "GABAergic Neurons 1", "GABAergic Progenitors 1")] = FALSE
Mito210_ARID1b_rep2$pseudoSub[Mito210_ARID1b_rep2$CellType %in% c("Newborn PNs", "IPC", "Newborn DL PNs")] = TRUE

saveRDS(Mito210_ARID1b_rep2, file = "celltypesSeur.rds")
saveRDS(Mito210_ARID1b_rep2, file = "clusteredSeur.rds")

setwd("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_1M_Harmony")
ARID1b_1M_harmony <- readRDS(file = "celltypesSeur.rds")

ARID1b_1M_harmony$pseudoSub = NA
ARID1b_1M_harmony$pseudoSub[ARID1b_1M_harmony$CellType %in% c("aRG", "Cycling Progenitors", "Subcortical",
                                                                  "Cortical Hem", "Cycling GABAergic Progenitors 1",
                                                                  "GABAergic Neurons 1", "GABAergic Progenitors 1")] = FALSE
ARID1b_1M_harmony$pseudoSub[ARID1b_1M_harmony$CellType %in% c("Newborn PNs", "IPC", "Newborn DL PNs")] = TRUE

saveRDS(ARID1b_1M_harmony, file = "celltypesSeur.rds")
saveRDS(ARID1b_1M_harmony, file = "Mito210_ARID1b_1M_harmonizedObj.rds")




#DEseq for pseudotime subsets of cells
#I made a metadata column in each object called "psuedoSub"(sic) with T/F based on whether that cell is included in the trajectory of interest for Fig. 3
for (dataset in c("Mito210_PTEN_1M_rep1/",
                  "Mito210_PTEN_1M_rep2/",
                  "Mito210_ARID1b_1M_scRNAseq_rep1_3.15.21/",
                  "Mito210_ARID1b_rep2_6.9.21/")) {
  
  #Load Seurat Object
  seur = readRDS(paste0(dataset,"celltypesSeur.rds"))
  
  #Set Ident to cluster or celltypes, whatever groups you want to seperate before detecting DEGs in that group
  Idents(seur) = "pseudoSub"
  id = TRUE
  
  #Set condition to the metadata column you want DEGs between, and base to the "wildtype" or base level of that column
  condition = "treat"
  base = "wt"
  
  #set combineOn to the metadata column that contains the samples (i.e. different organoids)
  combineOn = "orig.ident"
  
  dir=paste0(dataset,"wtvMut_DEGs_pseudoSub")
  
  #This loop will run DE analysis for each cluster and save a .xlsx file for each!
  degs <- combineDE(seur, id=id, condition=condition, base=base, combineOn=combineOn)
  degs$gene = rownames(degs)
  if (length(degs)>0) {
    write_xlsx(degs, path=paste0(dir,".DEGs.xlsx"), format_headers = T)
  }
}

setwd("~/Documents/Harvard/Arlotta Lab")
#DEseq for pseudotime subsets of cells for harmonized datasets
#I made a metadata column in each object called "psuedoSub"(sic) with T/F based on whether that cell is included in the trajectory of interest for Fig. 3
for (dataset in c("Mito210_ARID1b_1M_Harmony/",
                  "Mito210_PTEN_1M_Harmony/")) {
  
  #Load Seurat Object
  seur = readRDS(paste0(dataset,"celltypesSeur.rds"))
  
  #Set Ident to cluster or celltypes, whatever groups you want to seperate before detecting DEGs in that group
  Idents(seur) = "pseudoSub"
  id = TRUE
  
  #Set condition to the metadata column you want DEGs between, and base to the "wildtype" or base level of that column
  condition = "treat"
  base = "wt"
  
  #set combineOn to the metadata column that contains the samples (i.e. different organoids)
  combineOn = "orig.ident"
  
  dir=paste0(dataset,"wtvMut_DEGs_pseudoSub")
  
  #This loop will run DE analysis for each cluster and save a .xlsx file for each!
  degs <- combineDE(seur, id=id, condition=condition, base=base, combineOn=combineOn)
  degs$gene = rownames(degs)
  if (length(degs)>0) {
    write_xlsx(degs, path=paste0(dir,".DEGs.xlsx"), format_headers = T)
  }
}


#Read in DEGs for each and find overlaps and GO terms
allGenesUp = list()
allGenesDown = list()
datasets =  c("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_1M_scRNAseq_rep1_3.15.21/",
              "~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_rep2_6.9.21/",
              "~/Documents/Harvard/Arlotta Lab/Mito210_PTEN_1M_rep1/",
              "~/Documents/Harvard/Arlotta Lab/Mito210_PTEN_1M_rep2/")
names = c("ARID1b_1M_rep1", "ARID1b_1M_rep2","PTEN_1M_rep1","PTEN_1M_rep2")
dat = data.frame()
for (d in 1:length(datasets)) {
  dir=paste0(datasets[[d]],"wtvMut_DEGs_pseudoSub")
  xls = paste0(dir,".DEGs.xlsx")
  if (file.exists(xls)) {
    res = read_xlsx(xls)
    res = res[!is.na(res$pvalue),]
    res$dataset=names[[d]]
    genesUp = res[res$padj<0.05 & res$log2FoldChange>0,"gene"]
    genesDown = res[res$padj<0.05 & res$log2FoldChange<0,"gene"]
    allGenesUp[[names[[d]]]] = genesUp$gene
    allGenesDown[[names[[d]]]] = genesDown$gene
    dat = rbind(dat,res)
  }
}

ckU = compareCluster(allGenesUp, fun="enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP", pvalueCutoff=0.15)
ckD = compareCluster(allGenesDown, fun="enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP", pvalueCutoff=0.15)

ckU2 = clusterProfiler::simplify(ckU,cutoff=0.7, by="p.adjust", select_fun=min)
ckD2 = clusterProfiler::simplify(ckD,cutoff=0.7, by="p.adjust", select_fun=min)

dotplot(ckU2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  ggtitle("Upregulated in mutant") + scale_color_viridis(option="D", end=0.8)
ggsave("pseudo-GO.compareDatasets-Up.pdf")
write.csv(ckU2, file = "ckU_GOterms.csv", row.names = FALSE)

dotplot(ckD2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  ggtitle("Downregulated in mutant") + scale_color_viridis(end=0.8)
write.csv(ckD2, file = "ckD_GOterms.csv", row.names = FALSE)

ggsave("pseudo-GO.compareDatasets-Down.pdf")

#Harmonized objects
allGenesUp.Harmony = list()
allGenesDown.Harmony = list()
datasets =  c("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_1M_Harmony/",
              "~/Documents/Harvard/Arlotta Lab/Mito210_PTEN_1M_Harmony/")
names = c("ARID1b_1M", "PTEN_1M")
dat = data.frame()
for (d in 1:length(datasets)) {
  dir=paste0(datasets[[d]],"wtvMut_DEGs_pseudoSub")
  xls = paste0(dir,".DEGs.xlsx")
  if (file.exists(xls)) {
    res = read_xlsx(xls)
    res = res[!is.na(res$pvalue),]
    res$dataset=names[[d]]
    genesUp = res[res$padj<0.05 & res$log2FoldChange>0,"gene"]
    genesDown = res[res$padj<0.05 & res$log2FoldChange<0,"gene"]
    allGenesUp.Harmony[[names[[d]]]] = genesUp$gene
    allGenesDown.Harmony[[names[[d]]]] = genesDown$gene
    dat = rbind(dat,res)
  }
}

ckU.Harmony = compareCluster(allGenesUp.Harmony, fun="enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP", pvalueCutoff=0.15)
ckD.Harmony = compareCluster(allGenesDown.Harmony, fun="enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP", pvalueCutoff=0.15)

ckU2.Harmony = clusterProfiler::simplify(ckU.Harmony,cutoff=0.7, by="p.adjust", select_fun=min)
ckD2.Harmony = clusterProfiler::simplify(ckD.Harmony,cutoff=0.7, by="p.adjust", select_fun=min)

dotplot(ckU2.Harmony) + theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  ggtitle("Upregulated in mutant") + scale_color_viridis(option="D", end=0.8)
ggsave("pseudo-GO.compareDatasets-Up.pdf")
write.csv(ckU2.Harmony, file = "ckU_Harmony_GOterms.csv", row.names = FALSE)

dotplot(ckD2.Harmony) + theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  ggtitle("Downregulated in mutant") + scale_color_viridis(end=0.8)
write.csv(ckD2.Harmony, file = "ckD_Harmony_GOterms.csv", row.names = FALSE)

ggsave("pseudo-GO.compareDatasets-Down.pdf")


overlapGOtermsUp <- intersect()


#Get genes that overlap

allGenesUp.ARID1b <- c(allGenesUp$ARID1b_1M_rep1, allGenesUp$ARID1b_1M_rep2)
allGenesUp.PTEN <- c(allGenesUp$PTEN_1M_rep1, allGenesUp$PTEN_1M_rep2)

allGenesDown.ARID1b <- c(allGenesDown$ARID1b_1M_rep1, allGenesDown$ARID1b_1M_rep2)
allGenesDown.PTEN <- c(allGenesDown$PTEN_1M_rep1, allGenesDown$PTEN_1M_rep2)

overlap = intersect(allGenesUp.ARID1b, allGenesUp.PTEN)
overlap = union(overlap, intersect(allGenesDown.ARID1b, allGenesDown.PTEN))
datSlim = data.frame()
for (gene in overlap) {
  dats = dat[dat$gene==gene,]
  datSlim = rbind(datSlim, dats)
}

datSlim$dataset = factor(datSlim$dataset, levels=c("ARID1b_1M","PTEN_1M"))
dss = datSlim[datSlim$dataset=="ARID1b_1M",]
geneOrder = unique(dss[order(dss$log2FoldChange),]$gene)
datSlim$gene = factor(datSlim$gene, levels=geneOrder)
ggplot(datSlim, aes(x=gene, y=log2FoldChange, fill=log2FoldChange)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(dataset ~ ., scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)
ggsave("DEGs.ARID1BandPTENoverlap-rb.pdf")

overlap2 = intersect(allGenesUp$ARID1b_1M_rep1, allGenesUp$ARID1b_1M_rep2, allGenesUp$PTEN_1M_rep1, allGenesUp$PTEN_1M_rep2)
overlap2 = union(overlap2, intersect(allGenesDown$ARID1b_1M_rep1, allGenesDown$ARID1b_1M_rep2, allGenesDown$PTEN_1M_rep1, allGenesDown$PTEN_1M_rep2))
datSlim = data.frame()
for (gene in overlap2) {
  dats = dat[dat$gene==gene,]
  datSlim = rbind(datSlim, dats)
}
datSlim$dataset = factor(datSlim$dataset, levels=c("ARID1b_1M_rep1","PTEN_1M_rep1"))
dss = datSlim[datSlim$dataset=="ARID1b_1M_rep1",]
geneOrder = unique(dss[order(dss$log2FoldChange),]$gene)
datSlim$gene = factor(datSlim$gene, levels=geneOrder)
ggplot(datSlim, aes(x=gene, y=log2FoldChange, fill=log2FoldChange)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(dataset ~ ., scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)
ggsave("DEGs.ARID1BandPTENoverlap-rb.pdf")



overlap.Harmony = intersect(allGenesUp.Harmony$ARID1b_1M, allGenesUp.Harmony$PTEN_1M)
overlap.Harmony = union(overlap.Harmony, intersect(allGenesDown.Harmony$ARID1b_1M, allGenesDown.Harmony$PTEN_1M))
datSlim = data.frame()
for (gene in overlap.Harmony) {
  dats = dat[dat$gene==gene,]
  datSlim = rbind(datSlim, dats)
}
datSlim$dataset = factor(datSlim$dataset, levels=c("ARID1b_1M","PTEN_1M"))
dss = datSlim[datSlim$dataset=="ARID1b_1M",]
geneOrder = unique(dss[order(dss$log2FoldChange),]$gene)
datSlim$gene = factor(datSlim$gene, levels=geneOrder)
ggplot(datSlim, aes(x=gene, y=log2FoldChange, fill=log2FoldChange)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(dataset ~ ., scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)
ggsave("Harmony.DEGs.ARID1BandPTENoverlap-rb.pdf")




datSlim = data.frame()
for (gene in allGenesUp.Harmony$ARID1b_1M) {
  dats = dat[dat$gene==gene,]
  datSlim = rbind(datSlim, dats)
}
datSlim$dataset = factor(datSlim$dataset, levels=c("ARID1b_1M"))
dss = datSlim[datSlim$dataset=="ARID1b_1M",]
geneOrder = unique(dss[order(dss$log2FoldChange),]$gene)
datSlim$gene = factor(datSlim$gene, levels=geneOrder)
ggplot(datSlim, aes(x=gene, y=log2FoldChange, fill=log2FoldChange)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(dataset ~ ., scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)




