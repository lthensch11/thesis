# Load packages needed
library(dplyr)
library(Seurat)
library(ggplot2)
library(openxlsx)

setwd("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_rep2_6.9.21")

sessionInfo()

#Load in the proper scRNAseq file
Mito210_ARID1b_rep2 <- readRDS(file = "clusteredSeur.rds")
saveRDS(Mito210_ARID1b_rep2, file = "clusteredSeur.rds")

Mito210_ARID1b_rep2 <- readRDS(file = "celltypesSeur.rds")

#Visualize data with tSNE graphs/plots
DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", label = TRUE)
DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", label = TRUE, split.by = "treat")
DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", label = TRUE, split.by = "orig.ident", group.by = "treat")
DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", label = TRUE, split.by = "orig.ident", group.by = "orig.ident")
DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", label = TRUE, split.by = "orig.ident", group.by = "clusts_27PCs")
DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", label = TRUE, split.by = "orig.ident", group.by = "WTvMT.cluster.ids")

UMAP_Mito210_ARID1b_rep2 <- RunUMAP(Mito210_ARID1b_rep2, reduction = "pca", dims = 1:20)
DimPlot(UMAP_Mito210_ARID1b_rep2, reduction = "umap", label = TRUE, group.by = "clusts_27PCs")

FeaturePlot(UMAP_Mito210_ARID1b_rep2, reduction = "umap", features = c("DLX2", "DLX1", "GAD2"))

###Clusters 1, 16, 8, 14, 15, 20 more abundant in mut across replicates; 
###5, 2, 0, somewhat up in WT 

#WTvMT cluster assignment to generated split.by condition-specific tSNE plots
colnames(Mito210_ARID1b_rep2@meta.data)
Idents(Mito210_ARID1b_rep2) <- "clusts_27PCs" ## desired resolution needs to be active identity ie the resolution that you have used to assign cell types
Idents(Mito210_ARID1b_rep2) <- "orig.ident"
WTvMT.cluster.ids <- c("wt", "wt", "wt", "mut", "mut", "mut") #create a vector with new cell types in the same order ie cluster 1, position 1 # within c a fvector of names, a name for each cluster. for eg "aRG", "IP", .....
names(WTvMT.cluster.ids) <- levels(Mito210_ARID1b_rep2) #levels of the obj: identities that it has as active identities # with this function you assign the names to the clusters using the vector generated before
Mito210_ARID1b_rep2 <- RenameIdents(Mito210_ARID1b_rep2, WTvMT.cluster.ids) # for each cell, changing the name of the level assign to the cells to the new names
Mito210_ARID1b_rep2$WTvMT.cluster.ids <- Idents(Mito210_ARID1b_rep2)

#Find DEG's in the clusters
clustered_markers <- FindAllMarkers(Mito210_ARID1b_rep2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, max.cells.per.ident = 500)
markers <- FindAllMarkers(Mito210_ARID1b_rep2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,max.cells.per.ident = 500)

diffmark_WTvMut <- FindMarkers(Mito210_ARID1b_rep2, ident.1 = 'wt', ident.2 = 'mut', only.pos = TRUE, max.cells.per.ident = 500, group.by = "treat")
diffmark_WTvMut$gene <- rownames(diffmark_WTvMut)

?FindMarkers
DLX_neuron_markers <- FindMarkers(Mito210_ARID1b_rep2, ident.1 = "DLX neuron", only.pos = TRUE, max.cells.per.ident = 1000)
DLX_neuron_markers$gene <- rownames(DLX_neuron_markers)
DLX_aRG_markers <- FindMarkers(Mito210_ARID1b_rep2, ident.1 = "DLX aRG", only.pos = TRUE, max.cells.per.ident = 1000)
DLX_aRG_markers$gene <- rownames(DLX_aRG_markers)
DLX_aRGvneuron_markers <- FindMarkers(Mito210_ARID1b_rep2, ident.1 = "DLX neuron", ident.2 = "DLX aRG", only.pos = TRUE, max.cells.per.ident = 1000)
DLX_aRGvneuron_markers$gene <- rownames(DLX_aRGvneuron_markers)

##DEG's in WT and Mut conditions
wt_markers <- FindMarkers(Mito210_ARID1b_rep2, ident.1 = 'wt', onnly.pos = TRUE, max.cells.per.ident = 500, group.by = "treat")
wt_markers$gene <- rownames(wt_markers)

mut_markers <- FindMarkers(Mito210_ARID1b_rep2, ident.1 = 'mut', only.pos = TRUE, max.cells.per.ident = 500, group.by = "treat")
mut_markers$gene <- rownames(mut_markers)


#Graphical presentation of expressional profile
FeaturePlot(Mito210_ARID1b_rep2, features = c("DLX2", "DLX1", "GAD2"), reduction = "tsne")
DotPlot(Mito210_ARID1b_rep2, features = c("DLX2", "DLX1", "GAD2"))
VlnPlot(Mito210_ARID1b_rep2, features = c("DLX2", "DLX1", "GAD2"), pt.size = 0, ncol = 3)

FeaturePlot(Mito210_ARID1b_rep2, features = IN[, "V1"], reduction = "tsne")
#Clusters 14&16 express DLX2, 14 also express GAD2/DLX1

FeaturePlot(Mito210_ARID1b_rep2, features = Thalamus[, "V1"], reduction = "tsne")

FeaturePlot(Mito210_ARID1b_rep2, features = oRG[, "V1"], reduction = "tsne")

FeaturePlot(Mito210_ARID1b_rep2, features = Cycling[, "V1"], reduction = "tsne")

FeaturePlot(Mito210_ARID1b_rep2, features = IPCs[, "V1"], reduction = "tsne") 
#3,10,~9,12

FeaturePlot(Mito210_ARID1b_rep2, features = aRG[, "V1"], reduction = "tsne") 
#~0,1,2,4,11; not sure about 15/16

FeaturePlot(Mito210_ARID1b_rep2, features = Cortical.Hem[,"V1"], reduction = "tsne")

differentiation_markers <- c("FOXG1", "SNAP25", "BCL11B", "SATB2", "GAD1", "DLX2", "VIM", "PAX6", "SPARCL1", "EMX1", "NKX2-1", "EOMES", "TOP2A")
FeaturePlot(Mito210_ARID1b_rep2, features = differentiation_markers, reduction = "tsne")
FeaturePlot(Mito210_ARID1b_rep2, features = "FOXG1", reduction = "tsne")
DotPlot(Mito210_ARID1b_rep2, features = differentiation_markers)

FeaturePlot(Mito210_ARID1b_rep2, features = Excitatory.Neurons[, "V1"], reduction = "tsne") 
#7, ~3/6/9/13

FeaturePlot(Mito210_ARID1b_rep2, features = INM[, "V1"], reduction = "tsne") 
#12, ~4

cycling_genes <- c("MKI67", "TACC3", "TOP2A", "CDK1", "BIRC5", "CENPF")
FeaturePlot(Mito210_ARID1b_rep2, features = cycling_genes, reduction = "tsne")

cortical_hem_genes <- c("OTX2", "LMX1A", "WNT2B", "WNT5A", "WNT3A", "WLS", "FOXG1")
FeaturePlot(Mito210_ARID1b_rep2, features = cortical_hem_genes, reduction = "tsne")

IPC_genes <- c("EOMES", "PPP1R17", "PENK", "ELAVL4", "HES6", "NEUROD4") 
FeaturePlot(Mito210_ARID1b_rep2, features = IPC_genes, reduction = "tsne")

immature_neuronal_genes <- c("DCX", "NCAM1", "NEUROD1")
FeaturePlot(Mito210_ARID1b_rep2, features = immature_neuronal_genes, reduction = "tsne")

mature_neuronal_genes <- c("ENO2", "RBFOX3", "MAP2", "TUBB3", "NEFL", "NEFM", "NEFH", "GAP43")
FeaturePlot(Mito210_ARID1b_rep2, features = mature_neuronal_genes, reduction = "tsne")

synaptic_neuronal_genes <- c("DLG4", "SYP", "BSN")
FeaturePlot(Mito210_ARID1b_rep2, features = synaptic_neuronal_genes, reduction = "tsne")

FeaturePlot(Mito210_ARID1b_rep2, features = CFuPNs[,"V1"], reduction = "tsne")

FeaturePlot(Mito210_ARID1b_rep2, features = CPNs[,"V1"], reduction = "tsne")

FeaturePlot(Mito210_ARID1b_rep2, features = OLG[,"V1"], reduction = "tsne")

FeaturePlot(Mito210_ARID1b_rep2, features = "RELN", reduction = "tsne")


####CLUSTER IDENTIFICAITON SCRIPT####
#Data frame for relative prevalance of genes-related to specific cellular subtypes
Mito210_ARID1b_rep2_genetic_marker_rarity <- data.frame()

#anti-hem
antihem.table <- data.frame()

for (i in 1:lengths(anti.hem)[1]) {
  for(k in 0:20){
    antihem.table[anti.hem[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(anti.hem)[1]) {
  matches <- markers[markers$gene==anti.hem[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    antihem.table[anti.hem[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.antihem <- lengths(anti.hem)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["anti.hem", toString(i)] <- sum(antihem.table[,toString(i)]) / num.genes.antihem
}


# aRG
aRG.table <- data.frame()

for (i in 1:lengths(aRG)[1]) {
  for(k in 0:20){
    aRG.table[aRG[i, 1], toString(k)] <- 0
  }
}

# For each gene, find matching clusters
for (i in 1:lengths(aRG)[1]) {
  matches <- markers[markers$gene==aRG[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    aRG.table[aRG[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.aRG <- lengths(aRG)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["aRG", toString(i)] <- sum(aRG.table[,toString(i)]) / num.genes.aRG
}

#Astrocytes
astrocytes.table <- data.frame()

for (i in 1:lengths(Astrocytes)[1]) {
  for(k in 0:20){
    astrocytes.table[Astrocytes[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Astrocytes)[1]) {
  matches <- markers[markers$gene==Astrocytes[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    astrocytes.table[Astrocytes[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.astrocytes <- lengths(Astrocytes)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Astrocytes", toString(i)] <- sum(astrocytes.table[,toString(i)]) / num.genes.astrocytes
}

#CFuPN_CPNs aspecific
CFuPN.CPN.table <- data.frame()

for (i in 1:lengths(CFUPN_CPN.aspecific)[1]) {
  for(k in 0:20){
    CFuPN.CPN.table[CFUPN_CPN.aspecific[i, 1], toString(k)] <- 0
  }
}

# For each gene, find matching clusters
for (i in 1:lengths(CFUPN_CPN.aspecific)[1]) {
  matches <- markers[markers$gene==CFUPN_CPN.aspecific[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    CFuPN.CPN.table[CFUPN_CPN.aspecific[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.CFuPN.CPN <- lengths(CFUPN_CPN.aspecific)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["CFuPN.CPN", toString(i)] <- sum(CFuPN.CPN.table[,toString(i)]) / num.genes.CFuPN.CPN
}

#CFuPNs
CFuPN.table <- data.frame()

for (i in 1:lengths(CFuPNs)[1]) {
  for(k in 0:20){
    CFuPN.table[CFuPNs[i, 1], toString(k)] <- 0
  }
}

# For each gene, find matching clusters
for (i in 1:lengths(CFuPNs)[1]) {
  matches <- markers[markers$gene==CFuPNs[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    CFuPN.table[CFuPNs[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.CFuPNs <- lengths(CFuPNs)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["CFuPNs", toString(i)] <- sum(CFuPN.table[,toString(i)]) / num.genes.CFuPNs
}

#Cortical Hem
Cortical.hem.table <- data.frame()

for (i in 1:lengths(Cortical.Hem)[1]) {
  for(k in 0:20){
    Cortical.hem.table[Cortical.Hem[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Cortical.Hem)[1]) {
  matches <- markers[markers$gene==Cortical.Hem[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    Cortical.hem.table[Cortical.Hem[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.cortical.hem <- lengths(Cortical.Hem)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Cortical Hem", toString(i)] <- sum(CFuPN.CPN.table[,toString(i)]) / num.genes.cortical.hem
}


##Cycling progenitors
cycling.table <- data.frame()

for (i in 1:lengths(Cycling)[1]) {
  for(k in 0:20){
    cycling.table[Cycling[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Cycling)[1]) {
  matches <- markers[markers$gene==Cycling[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    cycling.table[Cycling[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.cycling <- lengths(Cycling)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Cycling Progenitors", toString(i)] <- sum(cycling.table[,toString(i)]) / num.genes.cycling
}

#Ectoderm
ectoderm.table <- data.frame()

for (i in 1:lengths(Ectoderm)[1]) {
  for(k in 0:20){
    ectoderm.table[Ectoderm[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Ectoderm)[1]) {
  matches <- markers[markers$gene==Ectoderm[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    ectoderm.table[Ectoderm[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.ectoderm <- lengths(Ectoderm)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Ectoderm", toString(i)] <- sum(ectoderm.table[,toString(i)]) / num.genes.ectoderm
}


##Excitatory neuronal progenitors
excitatory.neurons.table <- data.frame()

for (i in 1:lengths(Excitatory.Neurons)[1]) {
  for(k in 0:20){
    excitatory.neurons.table[Excitatory.Neurons[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Excitatory.Neurons)[1]) {
  matches <- markers[markers$gene==Excitatory.Neurons[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){
    excitatory.neurons.table[Excitatory.Neurons[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.excitatory.neurons <- lengths(Excitatory.Neurons)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Excitatory Neurons", toString(i)] <- sum(excitatory.neurons.table[,toString(i)]) / num.genes.excitatory.neurons
}

#Hypothalamus
hypothalamus.table <- data.frame()

for (i in 1:lengths(Hypothalamus)[1]) {
  for(k in 0:20){
    hypothalamus.table[Hypothalamus[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Hypothalamus)[1]) {
  matches <- markers[markers$gene==Hypothalamus[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    hypothalamus.table[Hypothalamus[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.hypothalamus <- lengths(Hypothalamus)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Hypothalamus", toString(i)] <- sum(hypothalamus.table[,toString(i)]) / num.genes.hypothalamus
}

#INM
INM.table <- data.frame()

for (i in 1:lengths(INM)[1]) {
  for(k in 0:20){
    INM.table[INM[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(INM)[1]) {
  matches <- markers[markers$gene==INM[i, 1],]
  if(lengths(matches)[1] <= 0){  
    next
  }
  for(k in 1:lengths(matches)[1]){  
    INM.table[INM[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.INM <- lengths(INM)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["INM", toString(i)] <- sum(INM.table[,toString(i)]) / num.genes.INM
}

#IPC
IPC.table <- data.frame()

for (i in 1:lengths(IPCs)[1]) {
  for(k in 0:20){
    IPC.table[IPCs[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(IPCs)[1]) {
  matches <- markers[markers$gene==IPCs[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    IPC.table[IPCs[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.IPC <- lengths(IPCs)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["IPCs", toString(i)] <- sum(IPC.table[,toString(i)]) / num.genes.IPC
}

#Melanocytes
melanocyte.table <- data.frame()

for (i in 1:lengths(Melanocytes)[1]) {
  for(k in 0:20){
    melanocyte.table[Melanocytes[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Melanocytes)[1]) {
  matches <- markers[markers$gene==Melanocytes[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    melanocyte.table[Melanocytes[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.melanocytes <- lengths(Melanocytes)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Melanocytes", toString(i)] <- sum(melanocyte.table[,toString(i)]) / num.genes.melanocytes
}

#Metabolic genes
metabolic.table <- data.frame()

for (i in 1:lengths(Metabolic.genes)[1]) {
  for(k in 0:20){
    metabolic.table[Metabolic.genes[i, 1], toString(k)] <- 0
  }
}

# For each gene, find matching clusters
for (i in 1:lengths(Metabolic.genes)[1]) {
  matches <- markers[markers$gene==Metabolic.genes[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    metabolic.table[Metabolic.genes[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.metabolic <- lengths(Metabolic.genes)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Metabolic Genes", toString(i)] <- sum(metabolic.table[,toString(i)]) / num.genes.metabolic
}

#Midbrain
midbrain.table <- data.frame()

for (i in 1:lengths(Midbrain)[1]) {
  for(k in 0:20){
    midbrain.table[Midbrain[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Midbrain)[1]) {
  matches <- markers[markers$gene==Midbrain[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    midbrain.table[Midbrain[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.midbrain <- lengths(Midbrain)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Midbrain", toString(i)] <- sum(midbrain.table[,toString(i)]) / num.genes.midbrain
}

#Neural Crest
neural.crest.table <- data.frame()

for (i in 1:lengths(Neural.crest)[1]) {
  for(k in 0:20){
    neural.crest.table[Neural.crest[i, 1], toString(k)] <- 0
  }
}

# For each gene, find matching clusters
for (i in 1:lengths(Neural.crest)[1]) {
  matches <- markers[markers$gene==Neural.crest[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    neural.crest.table[Neural.crest[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.neural.crest <- lengths(Neural.crest)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Neural Crest", toString(i)] <- sum(neural.crest.table[,toString(i)]) / num.genes.neural.crest
}

#Neuroepithelium
neuroepithelial.table <- data.frame()

for (i in 1:lengths(Neuroepithelial)[1]) {
  for(k in 0:20){
    neuroepithelial.table[Neuroepithelial[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Neuroepithelial)[1]) {
  matches <- markers[markers$gene==Neuroepithelial[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    neuroepithelial.table[Neuroepithelial[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.neuroeptihelial <- lengths(Neuroepithelial)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Neuroepithelial", toString(i)] <- sum(neuroepithelial.table[,toString(i)]) / num.genes.neuroeptihelial
}

#OLG
OLG.table <- data.frame()

for (i in 1:lengths(OLG)[1]) {
  for(k in 0:20){
    OLG.table[OLG[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(OLG)[1]) {
  matches <- markers[markers$gene==OLG[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    OLG.table[OLG[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.OLG <- lengths(OLG)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["OLG", toString(i)] <- sum(OLG.table[,toString(i)]) / num.genes.OLG
}


#oRG
oRG.table <- data.frame()

for (i in 1:lengths(oRG)[1]) {
  for(k in 0:20){
    oRG.table[oRG[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(oRG)[1]) {
  matches <- markers[markers$gene==oRG[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    oRG.table[oRG[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.oRG <- lengths(oRG)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["oRG", toString(i)] <- sum(oRG.table[,toString(i)]) / num.genes.oRG
}

#Retina
retina.table <- data.frame()

for (i in 1:lengths(Retina)[1]) {
  for(k in 0:20){
    retina.table[Retina[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Retina)[1]) {
  matches <- markers[markers$gene==Retina[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    retina.table[Retina[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.retina <- lengths(Retina)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Retina", toString(i)] <- sum(retina.table[,toString(i)]) / num.genes.retina
}

#Striatum
striatum.table <- data.frame()

for (i in 1:lengths(Striatum)[1]) {
  for(k in 0:20){
    striatum.table[Striatum[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Striatum)[1]) {
  matches <- markers[markers$gene==Striatum[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    striatum.table[Striatum[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.striatum <- lengths(Striatum)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Striatum", toString(i)] <- sum(striatum.table[,toString(i)]) / num.genes.striatum
}

#Subplate
subplate.table <- data.frame()

for (i in 1:lengths(Subplate)[1]) {
  for(k in 0:20){
    subplate.table[Subplate[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Subplate)[1]) {
  matches <- markers[markers$gene==Subplate[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    subplate.table[Subplate[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.subplate <- lengths(Subplate)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Subplate", toString(i)] <- sum(subplate.table[,toString(i)]) / num.genes.subplate
}

#Synpatic
synaptic.table <- data.frame()

for (i in 1:lengths(Synaptic)[1]) {
  for(k in 0:20){
    synaptic.table[Synaptic[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Synaptic)[1]) {
  matches <- markers[markers$gene==Synaptic[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    synaptic.table[Synaptic[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.synaptic <- lengths(Synaptic)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Synaptic", toString(i)] <- sum(synaptic.table[,toString(i)]) / num.genes.synaptic
}

#Thalamus
thalamus.table <- data.frame()

for (i in 1:lengths(Thalamus)[1]) {
  for(k in 0:20){
    thalamus.table[Thalamus[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Thalamus)[1]) {
  matches <- markers[markers$gene==Thalamus[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    thalamus.table[Thalamus[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.thalamus <- lengths(Thalamus)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Thalamus", toString(i)] <- sum(thalamus.table[,toString(i)]) / num.genes.thalamus
}

#Ventral mix
ventral.mix.table <- data.frame()

for (i in 1:lengths(Ventral.mix)[1]) {
  for(k in 0:20){
    ventral.mix.table[Ventral.mix[i, 1], toString(k)] <- 0
  }
}

for (i in 1:lengths(Ventral.mix)[1]) {
  matches <- markers[markers$gene==Ventral.mix[i, 1],]
  if(lengths(matches)[1] <= 0){  # If no matches, move on
    next
  }
  for(k in 1:lengths(matches)[1]){  # If match, then add 1 to that column's cell
    ventral.mix.table[Ventral.mix[i, 1], toString(matches[k, "cluster"])] = 1
  }
}

num.genes.ventral.mix <- lengths(Ventral.mix)[1]
for(i in 0:20){
  Mito210_ARID1b_rep2_genetic_marker_rarity["Ventral Mix", toString(i)] <- sum(ventral.mix.table[,toString(i)]) / num.genes.ventral.mix
}


write.xlsx(Mito210_ARID1b_rep2_genetic_marker_rarity, "~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_rep2_6.9.21/Mito210_ARID1b_rep2_1M_cluster_ident_frequency.xlsx", row.names = TRUE)


#####Assingment of new cluster IDs based on rarity output####
Idents(Mito210_ARID1b_rep2) <- "clusts_27PCs"
new.clusters.id <- c("aRG", "aRG", "Cycling aRG", "Excitatory neurons", "Cycling aRG", "Excitatory neurons", 
                     "Excitatory neurons", "Excitatory neurons", "DLX neuron", "Excitatory neurons", "IPC",
                     "Cycling aRG", "Cycling IPC", "Excitatory neurons", "DLX neuron", "Cycling aRG", 
                     "DLX aRG", "Cortical Hem", "Excitatory neurons", "Cycling IPC", "DLX neuron")
names(new.clusters.id) <- levels(Mito210_ARID1b_rep2)
Mito210_ARID1b_rep2 <- RenameIdents(Mito210_ARID1b_rep2, new.clusters.id)

DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", label = TRUE, pt.size = 0.5)
DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "orig.ident")
DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "treat")

clusters.id.v2 <- c("aRG","aRG","Cycling Progenitors","Newborn PNs","Cycling Progenitors",
                    "Newborn DL PNs","Newborn PNs","Newborn DL PNs","GABAergic Neurons 1",
                    "Newborn PNs","IPC","Cycling Progenitors","Cycling Progenitors","Newborn DL PNs",
                    "GABAergic Neurons 1","Cycling GABAergic Progenitors 1","GABAergic Progenitors 1",
                    "Cortical Hem","Subcortical","Cycling Progenitors","GABAergic Neurons 1")
names(clusters.id.v2) <- levels(Mito210_ARID1b_rep2)
Mito210_ARID1b_rep2 <- RenameIdents(Mito210_ARID1b_rep2, clusters.id.v2)

#Create a new CellType meta-data column
Mito210_ARID1b_rep2$CellType = NA
Mito210_ARID1b_rep2$CellType[Mito210_ARID1b_rep2$clusts_27PCs %in% c(0,1)] = "aRG" 
Mito210_ARID1b_rep2$CellType[Mito210_ARID1b_rep2$clusts_27PCs %in% c(2,4,11,12,19)] = "Cycling Progenitors" 
Mito210_ARID1b_rep2$CellType[Mito210_ARID1b_rep2$clusts_27PCs %in% c(3,6,9)] = "Newborn PNs" 
Mito210_ARID1b_rep2$CellType[Mito210_ARID1b_rep2$clusts_27PCs %in% c(5,7,13)] = "Newborn DL PNs" 
Mito210_ARID1b_rep2$CellType[Mito210_ARID1b_rep2$clusts_27PCs %in% c(8,14,20)] = "GABAergic Neurons 1" 
Mito210_ARID1b_rep2$CellType[Mito210_ARID1b_rep2$clusts_27PCs %in% c(15)] = "Cycling GABAergic Progenitors 1" 
Mito210_ARID1b_rep2$CellType[Mito210_ARID1b_rep2$clusts_27PCs %in% c(16)] = "GABAergic Progenitors 1"
Mito210_ARID1b_rep2$CellType[Mito210_ARID1b_rep2$clusts_27PCs %in% c(10)] = "IPC" 
Mito210_ARID1b_rep2$CellType[Mito210_ARID1b_rep2$clusts_27PCs %in% c(17)] = "Cortical Hem" 
Mito210_ARID1b_rep2$CellType[Mito210_ARID1b_rep2$clusts_27PCs %in% c(18)] = "Subcortical" 


levels_A2 =  c("aRG","Cycling Progenitors","Cortical Hem","IPC",
            "Newborn DL PNs","Newborn PNs","Subcortical", 
            "GABAergic Progenitors 1","Cycling GABAergic Progenitors 1","GABAergic Neurons 1")
cols_A2 = c("#41ae76","#bbcc33", "#bebada", "#fdb462",
         "#f768a1","#fa9fb5","#77aadd",
         "turquoise","aquamarine","turquoise4")

colsA2 = cols_A2[match(levels(factor((Mito210_ARID1b_rep2)$CellType)),levels_A2)]

DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", cols=colsA2, group.by = "CellType") + ggtitle("")+ NoAxes() + NoLegend()
DimPlot(Mito210_ARID1b_rep2, reduction="tsne", cols=colsA2, group.by = "CellType", split.by = "treat") + NoAxes() + ggtitle("")

Idents(Mito210_ARID1b_rep2) <- "CellType"

#To align Martina's cluster identities with my own. Correction was needed to have the cell names match
#the one's from Martina's Seurat object.
a <- substr(rownames(Mito210_ARID1b_rep2@meta.data), 1, nchar(rownames(Mito210_ARID1b_rep2@meta.data))-2)
Mito210_ARID1b_rep2$MartinaCellTypes = Martina.ARID1b.rep2$CellType[rownames(Mito210_ARID1b_rep2@meta.data)]

p1 = DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", group.by = "MartinaCellTypes", label = TRUE)
p2 = DimPlot(Mito210_ARID1b_rep2, label = TRUE)

p1+p2



####GO Analysis####
#load all the libraries
library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(org.Hs.eg.db)
library(readr)
library(clusterProfiler)


#GO graph
symbols <- markers[,"gene"]
#c(row.names(markers_SUV35Excitatory))
#c(correlating_genes$Row.names)
symbols <- unlist(symbols)
ego = enrichGO(gene= symbols, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP", pvalueCutoff=0.6)
dotplot(ego) + theme(axis.text.x=element_text(angle=45, hjust=1)) +theme(plot.title = element_text(hjust = 0.5))+ ggtitle('GO Mito210 ARID1B rep2 6.9.21')

DEG_ARID1B_markers <- diffmark_WTvMut[, "gene"]
ego_DEG_ARID1b = enrichGO(gene = DEG_ARID1B_markers, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.6)
dotplot(ego_DEG_ARID1b) + theme(axis.text.x=element_text(angle=45, hjust=1)) +theme(plot.title = element_text(hjust = 0.5))+ ggtitle('GO ARID1B DEG')

cluster8_ARID1b_markers <- markers[6021:6422, "gene"]
ego_cluster8 = enrichGO(gene = cluster8_ARID1b_markers, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                        ont = "BP", pvalueCutoff = 0.6)
dotplot(ego_cluster8) + theme(axis.text.x=element_text(angle=45, hjust=1)) +theme(plot.title = element_text(hjust = 0.5))+ ggtitle('GO Mito210 ARID1B DLX Neuron Mut (8)')

cluster14_ARID1b_markers <- markers[9928:10190, "gene"]
head(cluster14_ARID1b_markers)
tail(cluster14_ARID1b_markers)
ego_cluster14 = enrichGO(gene = cluster14_ARID1b_markers, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.6)
dotplot(ego_cluster14) + theme(axis.text.x=element_text(angle=45, hjust=1)) +theme(plot.title = element_text(hjust = 0.5))+ ggtitle('GO Mito210 ARID1B DLX Neuron Mut (14)')


wt_DEG <- wt_markers[, "gene"]
ego_wt_DEG = enrichGO(gene = wt_DEG, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",ont = "BP", pvalueCutoff = 0.6)
dotplot(ego_wt_DEG) + theme(axis.text.x=element_text(angle=45, hjust=1)) +theme(plot.title = element_text(hjust = 0.5))+ ggtitle('GO Mito210 ARID1B WT DEG')

mut_DEG <- mut_markers[, "gene"]
ego_mut_DEG = enrichGO(gene = mut_DEG, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",ont = "BP", pvalueCutoff = 0.6)
dotplot(ego_mut_DEG) + theme(axis.text.x=element_text(angle=45, hjust=1)) +theme(plot.title = element_text(hjust = 0.5))+ ggtitle('GO Mito210 ARID1B Mut DEG')

DLX_neuron_genes <- DLX_neuron_markers[,"gene"]
ego_DLX_neuron = enrichGO(gene = DLX_neuron_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.6)
dotplot(ego_DLX_neuron) + theme(axis.text.x=element_text(angle=45, hjust=1)) +theme(plot.title = element_text(hjust = 0.5))+ ggtitle('GO Mito210 ARID1B DLX Neuron')

DLX_aRG_genes <- DLX_aRG_markers[,"gene"]
ego_DLX_aRG = enrichGO(gene = DLX_aRG_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.6)
dotplot(ego_DLX_aRG) + theme(axis.text.x=element_text(angle=45, hjust=1)) +theme(plot.title = element_text(hjust = 0.5))+ ggtitle('GO Mito210 ARID1B DLX aRG')



setwd("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_rep2_6.9.21/DESeq2 DEGs Cell Type/")
res = read_xlsx("Newborn DL PNs.DEGs.xlsx")
res = res[!is.na(res$pvalue),]
NewbornDLPN.genesUp = res[res$padj<0.05 & res$log2FoldChange>0,"gene"]
NewbornDLPN.genesUp <- unlist(NewbornDLPN.genesUp)
NewbornDLPN.genesDown = res[res$padj<0.05 & res$log2FoldChange<0,"gene"]
NewbornDLPN.genesDown <- unlist(NewbornDLPN.genesDown)

ego_NewbornDLPN_genesUp = enrichGO(gene = NewbornDLPN.genesUp, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
simplify_ego_NewbornDLPN_genesUp <- simplify(ego_NewbornDLPN_genesUp, cutoff = 0.7, by="p.adjust", select_fun=min)
dotplot(simplify_ego_NewbornDLPN_genesUp) + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle('GO ARID1b 1M rep1 Newborn DL PNs genesUp DEG')

ego_NewbornDLPN_genesDown = enrichGO(gene = NewbornDLPN.genesDown, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.25)
simplify_ego_NewbornDLPN_genesDown <- simplify(ego_NewbornDLPN_genesDown, cutoff = 0.7, by="p.adjust", select_fun=min)
dotplot(simplify_ego_NewbornDLPN_genesDown) + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle('GO ARID1b 1M rep1 Newborn DL PNs genesDown DEG')

write.csv(data.frame(simplify_ego_NewbornDLPN_genesUp), file = "NewbornDLPN_genesUp_DEGs_GOterms.csv", row.names = FALSE)
write.csv(data.frame(simplify_ego_NewbornDLPN_genesDown), file = "NewbornDLPN_genesDown_DEGs_GOterms.csv", row.names = FALSE)


####KEGG ANALYSIS####
#KEGG graph
#converting the list of genes in Entrez ID.the function requires a character vector starting from a matrix
#symbols <- row.names(markers_PTEN)
#symbols <- all_DEP_1[,1]
#symbols <- as.matrix(symbols)
#symbols <- as.vector(symbols)


#what input mapIDs requires??
?mapIds
# use mapIds method to obtain Entrez IDs
markers_mapIds <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
kk <-  enrichKEGG(gene= markers_mapIds, organism= "hsa", pvalueCutoff = 0.5)
dotplot(kk) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle('correlating genes\n no filter KEGG all Markers')


DEG_WTvMut_ARID1b_markers_mapIds <- mapIds(org.Hs.eg.db, DEG_ARID1B_markers, 'ENTREZID', 'SYMBOL')
kk_DEG_WTvMut <- enrichKEGG(gene = DEG_WTvMut_ARID1b_markers_mapIds, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_DEG_WTvMut) + theme(axis.text.x = element_text(angle=45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG DEG WT v. Mut')

wt_DEG_markers_mapIds <- mapIds(org.Hs.eg.db, wt_DEG, 'ENTREZID', 'SYMBOL')
kk_wt_DEG <- enrichKEGG(gene = wt_DEG_markers_mapIds, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_wt_DEG) + theme(axis.text.x = element_text(angle=45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG DEG WT')

mut_DEG_markers_mapIds <- mapIds(org.Hs.eg.db, mut_DEG, 'ENTREZID', 'SYMBOL')
kk_mut_DEG <- enrichKEGG(gene = mut_DEG_markers_mapIds, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_mut_DEG) + theme(axis.text.x = element_text(angle=45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG DEG Mut')

DLX_neuron_markers_mapIds <- mapIds(org.Hs.eg.db, DLX_neuron_genes, 'ENTREZID', 'SYMBOL')
kk_DLX_neuron <- enrichKEGG(gene = DLX_neuron_markers_mapIds, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_DLX_neuron) + theme(axis.text.x = element_text(angle=45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG DLX Neuron')

DLX_aRG_markers_mapIds <- mapIds(org.Hs.eg.db, DLX_aRG_genes, 'ENTREZID', 'SYMBOL')
kk_DLX_aRG <- enrichKEGG(gene = DLX_aRG_markers_mapIds, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_DLX_aRG) + theme(axis.text.x = element_text(angle=45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG DLX aRG')


#symbols <- as.matrix(symbols)
#symbols <- as.vector(symbols)


diffmark_aRG <- as.vector(diffmark_aRG)
diffmark_neuronal <- as.vector(diffmark_neuronal)



####Monocle3/Pseudotime Constructs####
#Load Monocle
library(monocle3)
library(Seurat)

#Helper function! This lets you pick the start point of your pseudotime. 
#Here I pick based on the "ScinaCellType" metadata column, and set the default start point to be "Cycling" cells - change as you like!
get_earliest_principal_node <- function(cds,cell_type="Cycling Progenitors"){
  cell_ids <- which(colData(cds)[, "CellType"] == cell_type)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

#Load data first
fd = data.frame("gene_short_name" = rownames(Mito210_ARID1b_rep2))
rownames(fd) = rownames(Mito210_ARID1b_rep2)
cds <- new_cell_data_set(GetAssayData(Mito210_ARID1b_rep2,slot="counts"),
                         cell_metadata = Mito210_ARID1b_rep2@meta.data,
                         gene_metadata = fd)

#Subset so control cells equals mutant cells
table(colData(cds)$treat) #look at which one has fewer cells and downsample the other to that number
wtCells = rownames(colData(cds)[colData(cds)$treat=="wt",])
mutCells = sample(rownames(colData(cds)[colData(cds)$treat=="mut",]), 16434,replace = F)
cds = cds[,c(wtCells,mutCells)]

#Run Monocle
cds <- preprocess_cds(cds, num_dim = 30)
#plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,umap.fast_sgd=TRUE) #makes the UMAP
cds <- cluster_cells(cds, partition_qval=0.5) #playing with "partition_qval" from 0-1 will change how many different lineages your cells are split into
cds <- learn_graph(cds)

levels =  c("aRG","Cajal-Retzius","Cortical Hem","Cycling Progenitors","IPC",
            "Newborn DL PNs","Newborn PNs","Subcortical", "Unknown", "Choroid Plexus/Cortical Hem",
            "CFuPNs","CPNs","oRG","PNs","GABAergic Neurons 2",
            "oRG/Astroglia","Astroglia","GABAergic Progenitors 2", "Cycling GABAergic Progenitors 2",
            "GABAergic Progenitors 1","Cycling GABAergic Progenitors 1","GABAergic Neurons 1")

cols = c("#41ae76","#ee8866","#bebada","#bbcc33","#fdb462",
         "#f768a1","#fa9fb5","#77aadd","darkgray", "powderblue",
         "#cc6677","#882255","#225522","#aa4499","#332288",
         "#009988","#0077bb","#5B65AE","#B87ACF",
         "turquoise","aquamarine","turquoise4")
cols = cols[match(levels(factor(colData(cds)$CellType)),levels)]

p1 = plot_cells(cds, color_cells_by = "CellType", 
                cell_size = 0.5,
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                labels_per_group = 0) + scale_color_manual(values=cols) + theme_void()
cds = order_cells(cds,root_pr_nodes=get_earliest_principal_node(cds)) #Gives each cell a pseudotime value. Uses the helper function from above to set start point

colData(cds)$clusts_27PCs = factor(colData(cds)$clusts_27PCs)

#remove helper funnction here to be able to select starting node for trajectory analysis
#run but this time subset out by progenitors v neuron fates
p2=plot_cells(cds, color_cells_by = "clusts_27PCs",show_trajectory_graph = F,
              label_groups_by_cluster = F,group_label_size = 4,cell_size = 0.5) + theme_void()
p3 = plot_cells(cds,show_trajectory_graph = T,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=1.5,
                cell_size = 0.5) + theme_void()
p4=plot_cells(cds, color_cells_by = "treat", show_trajectory_graph = F, label_cell_groups = F,cell_size = 0.5) + theme_void()

p1
p2
p3
p4



colData(cds)$pseudotime <- pseudotime(cds)
df = as.data.frame(colData(cds))
df = df[is.finite(df$pseudotime),]

ggplot(df, aes(x=pseudotime, y=CellType)) +
  geom_density_ridges(aes(fill=CellType)) +
  theme_classic() + scale_color_manual(values = cols) +
  ggtitle('ARID1b 1M rep2 Pseudotime Original Split.by CellType')

ggplot(df, aes(x=pseudotime, y=treat)) +
  geom_density_ridges(aes(fill=treat)) +
  theme_void() +
  ggtitle('ARID1b 1M rep2 Pseudotime Original Split.by Treatment')

#to subset(object, subset = metadata name %in% c(all cell types interested in)
#only cortical structures included#
pseudo_subset_ARID1b_rep2 <- subset(Mito210_ARID1b_rep2, CellType %in% c("aRG", "Cajal-Retzius", "Cycling GABAergic Progenitors 1", "GABAergic Neurons 1", "GABAergic Progenitors 1", "Cycling Progenitors", "IPC", "Newborn PNs", "Newborn DL PNs"))

#subset by treatment groups#
pseudo_WT_ARID1b_rep2 <- subset(Mito210_ARID1b_rep2, treat %in% c("wt"))
pseudo_MT_ARID1b_rep2 <- subset(Mito210_ARID1b_rep2, treat %in% c("mut"))


#WT#
fd_pseudo_WT_ARID1b_rep2 = data.frame("gene_short_name" = rownames(pseudo_WT_ARID1b_rep2))
rownames(fd_pseudo_WT_ARID1b_rep2) = rownames(pseudo_WT_ARID1b_rep2)
cds_pseudo_WT_ARID1b_rep2 <- new_cell_data_set(GetAssayData(pseudo_WT_ARID1b_rep2,slot="counts"),
                         cell_metadata = pseudo_WT_ARID1b_rep2@meta.data,
                         gene_metadata = fd_pseudo_WT_ARID1b_rep2)

#Run Monocle
cds_pseudo_WT_ARID1b_rep2 <- preprocess_cds(cds_pseudo_WT_ARID1b_rep2, num_dim = 30)
#plot_pc_variance_explained(cds)
cds_pseudo_WT_ARID1b_rep2 <- reduce_dimension(cds_pseudo_WT_ARID1b_rep2,umap.fast_sgd=TRUE) #makes the UMAP
p1=plot_cells(cds_pseudo_WT_ARID1b_rep2, color_cells_by = "CellType") #choose the metadata column you're interested in
cds_pseudo_WT_ARID1b_rep2 <- cluster_cells(cds_pseudo_WT_ARID1b_rep2, partition_qval=1) #playing with "partition_qval" from 0-1 will change how many different lineages your cells are split into
cds_pseudo_WT_ARID1b_rep2 <- learn_graph(cds_pseudo_WT_ARID1b_rep2)

#remove helper funnction here to be able to select starting node for trajectory analysis
#run but this time subset out by progenitors v neuron fates
cds_pseudo_WT_ARID1b_rep2 = order_cells(cds_pseudo_WT_ARID1b_rep2,root_pr_nodes=get_earliest_principal_node(cds_pseudo_WT_ARID1b_rep2)) #Gives each cell a pseudotime value. Uses the helper function from above to set start point
p2=plot_cells(cds_pseudo_WT_ARID1b_rep2,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=1.5)

p1+p2

colData(cds_pseudo_WT_ARID1b_rep2)$pseudotime <- pseudotime(cds_pseudo_WT_ARID1b_rep2)
df = as.data.frame(colData(cds_pseudo_WT_ARID1b_rep2))
df = df[is.finite(df$pseudotime),]

ggplot(df, aes(x=pseudotime, y=CellType)) +
  geom_density_ridges(aes(fill=CellType)) +
  theme_classic() +
  ggtitle('ARID1b 1M rep2 Pseudotime WT Split.by CellType')



#Mut#
fd_pseudo_MT_ARID1b_rep2 = data.frame("gene_short_name" = rownames(pseudo_MT_ARID1b_rep2))
rownames(fd_pseudo_MT_ARID1b_rep2) = rownames(Mito210_ARID1b_rep2)
cds_pseudo_MT_ARID1b_rep2 <- new_cell_data_set(GetAssayData(Mito210_ARID1b_rep2,slot="counts"),
                         cell_metadata = Mito210_ARID1b_rep2@meta.data,
                         gene_metadata = fd_pseudo_MT_ARID1b_rep2)


#Run Monocle
cds_pseudo_MT_ARID1b_rep2 <- preprocess_cds(cds_pseudo_MT_ARID1b_rep2, num_dim = 30)
#plot_pc_variance_explained(cds)
cds_pseudo_MT_ARID1b_rep2 <- reduce_dimension(cds_pseudo_MT_ARID1b_rep2,umap.fast_sgd=TRUE) #makes the UMAP
p1=plot_cells(cds_pseudo_MT_ARID1b_rep2, color_cells_by = "CellType") #choose the metadata column you're interested in
cds_pseudo_MT_ARID1b_rep2 <- cluster_cells(cds_pseudo_MT_ARID1b_rep2, partition_qval=1) #playing with "partition_qval" from 0-1 will change how many different lineages your cells are split into
cds_pseudo_MT_ARID1b_rep2 <- learn_graph(cds_pseudo_MT_ARID1b_rep2)

#remove helper funnction here to be able to select starting node for trajectory analysis
#run but this time subset out by progenitors v neuron fates
cds_pseudo_MT_ARID1b_rep2 = order_cells(cds_pseudo_MT_ARID1b_rep2,root_pr_nodes=get_earliest_principal_node(cds_pseudo_MT_ARID1b_rep2)) #Gives each cell a pseudotime value. Uses the helper function from above to set start point
p2=plot_cells(cds_pseudo_MT_ARID1b_rep2,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=1.5)

p1 + p2

colData(cds_pseudo_MT_ARID1b_rep2)$pseudotime <- pseudotime(cds_pseudo_MT_ARID1b_rep2)
df = as.data.frame(colData(cds_pseudo_MT_ARID1b_rep2))
df = df[is.finite(df$pseudotime),]

ggplot(df, aes(x=pseudotime, y=CellType)) +
  geom_density_ridges(aes(fill=CellType)) +
  theme_classic() +
  ggtitle('ARID1b 1M rep2 Pseudotime MT Split.by CellType')


DimPlot(Mito210_ARID1b_rep2, split.by = "treat")


####DESeq2 Analysis####
#Takes a Seurat Object and performs differential expression analysis 
#For each cluster/cell type
#By splitting the object into separate samples (i.e organoids) and using DESeq2,
#which accounts for noise between samples and then looks for DEGs that overcome that noise
install.packages("DESeq2")
library(Seurat)
library(DESeq2)
library(writexl)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rnaseqGene")

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
Mito210_ARID1b_rep2 = readRDS(file = "clusteredSeur.rds")

#Set Ident to cluster or celltypes, whatever groups you want to seperate before detecting DEGs in that group
Idents(Mito210_ARID1b_rep2) = "CellType"
id = "CellType"

#Set condition to the metadata column you want DEGs between, and base to the "wildtype" or base level of that column
condition = "treat"
base = "wt"

#set combineOn to the metadata column that contains the samples (i.e. different organoids)
combineOn = "orig.ident"

#This loop will run DE analysis for each cluster and save a .xlsx file for each!
for (id in levels(Mito210_ARID1b_rep2@active.ident)) {
  print(id)
  degs <- combineDE(Mito210_ARID1b_rep2, id=id, condition=condition, base=base, combineOn=combineOn)
  degs$gene = rownames(degs)
  if (length(degs)>0) {
    write_xlsx(degs, path=paste0(id,".DEGs.xlsx"), format_headers = T)
  }
}

##All GABAergic clusters were not considered because of their variable naming. 
##--> we may want to create an object with a combined GABAergic group

Idents(Mito210_ARID1b_rep2) = "clusts_28PCs"
id = "clusts_28PCs"
condition = "treat"
base = "wt"
combineOn = "orig.ident"

for (id in levels(Mito210_ARID1b_rep2@active.ident)) {
  print(id)
  degs <- combineDE(Mito210_ARID1b_rep2, id=id, condition=condition, base=base, combineOn=combineOn)
  degs$gene = rownames(degs)
  if (length(degs)>0) {
    write_xlsx(degs, path=paste0("clus",id,".DEGs.xlsx"), format_headers = T)
  }
}

##Since there was such a distinct wt/mt phenotype, many of the clusters did not get a deg list. 
##For example all of the GABAergic clusters had no noticeable cells in the WT orgs, so there was no deg generated to compare wt and mut groups


####Speckle/Cell Type Proportion####
library(devtools)
library(BiocManager)
BiocManager::install("org.Mm.eg.db")

devtools::install_github("Oshlack/speckle")
# devtools/remotes won't install Suggested packages from Bioconductor
BiocManager::install(c("CellBench", "BiocStyle", "scater"))

remotes::install_github("Oshlack/speckle", build_vignettes = TRUE, 
                        dependencies = "Suggest")

library(speckle)
library(limma)
library(ggplot2)


# Get some example data which has two groups, three cell types and two 
# biological replicates in each group
cell_info <- Mito210_ARID1b_rep2
head(cell_info)

# Run propeller testing for cell type proportion differences between the two 
# groups
?propeller
propeller.test.ARID1b.rep2 <- propeller(clusters = Mito210_ARID1b_rep2$CellType, sample = Mito210_ARID1b_rep2$orig.ident, 
                            group = Mito210_ARID1b_rep2$treat)
write.xlsx(propeller.test.ARID1b.rep2, file = "ARID1b_1M_rep2_propeller.xlsx")


# Plot cell type proportions
plotCellTypeProps(clusters=Mito210_ARID1b_rep2$CellType, sample=Mito210_ARID1b_rep2$orig.ident) + theme_classic() + theme(legend.position = "bottom") + scale_color_manual(values = colsA2)
plotCellTypeProps(clusters=Mito210_ARID1b_rep2$CellType, sample=Mito210_ARID1b_rep2$treat) + theme_classic()


a <- subset(Mito210_ARID1b_rep2, CellType %in% c("Cycling GABAergic Progenitors 1", "GABAergic Neurons 1", "GABAergic Progenitors 1", "Cycling Progenitors"))

plotCellTypeProps(clusters=a$CellType, sample=a$orig.ident, cols=cols_A2)

?plotCellTypeProps

b <- subset(Mito210_ARID1b_rep2, CellType %in% c("GABAergic Neurons 1"))
plotCellTypeProps(clusters=b$CellType, sample=b$orig.ident)




