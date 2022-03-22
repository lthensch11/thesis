# Load packages needed
library(dplyr)
library(Seurat)
library(ggplot2)
library(openxlsx)

install.packages("ggplot2", "openxlsx")

setwd("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_1M_scRNAseq_rep1_3.15.21")

sessionInfo()

#Load in the proper scRNAseq file
ARID1b_1M_rep1 <- readRDS(file = "clusteredSeur_ARID1B_1M.rds")
ARID1b_1M_rep1 <- readRDS(file = "celltypesSeur.rds")

#Visualize tSNE plot of ARID1b data
DimPlot(ARID1b_1M_rep1, reduction = "tsne", label = TRUE) + NoLegend()
DimPlot(ARID1b_1M_rep1, reduction = "tsne", label = TRUE, split.by = "orig.ident")
DimPlot(ARID1b_1M_rep1, reduction = "tsne", label = TRUE, split.by = "WTvMT.cluster.ids")
DimPlot(ARID1b_1M_rep1, reduction = "tsne", label = TRUE, group.by = "WTvMT.cluster.ids", split.by = "WTvMT.cluster.ids")

colnames(ARID1b_1M_rep1@meta.data)
Idents(ARID1b_1M_rep1) <- "clusts_28PCs" ## desired resolution needs to be active identity ie the resolution that you have used to assign cell types
Idents(ARID1b_1M_rep1) <- "orig.ident"
WTvMT.cluster.ids <- c("WT", "WT", "WT", "MT", "MT", "MT") #create a vector with new cell types in the same order ie cluster 1, position 1 # within c a fvector of names, a name for each cluster. for eg "aRG", "IP", .....
names(WTvMT.cluster.ids) <- levels(ARID1b_1M_rep1) #levels of the obj: identities that it has as active identities # with this function you assign the names to the clusters using the vector generated before
ARID1b_1M_rep1 <- RenameIdents(ARID1b_1M_rep1, WTvMT.cluster.ids) # for each cell, changing the name of the level assign to the cells to the new names
ARID1b_1M_rep1$WTvMT.cluster.ids <- Idents(ARID1b_1M_rep1)

markers <- FindAllMarkers(ARID1b_1M_rep1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,max.cells.per.ident = 500)

diffmark <- FindMarkers(ARID1b_1M_rep1, ident.1 = 1, ident.2 = 2, only.pos = TRUE, max.cells.per.ident = 500, group.by = "orig.ident")

cluster0 <- FindMarkers(ARID1b_1M_rep1, ident.1 = 0, only.pos = TRUE, max.cells.per.ident = 500, min.pct = 0.1)

diffmark_WTvMT <- FindMarkers(ARID1b_1M_rep1, ident.1 = c(1,2,3), ident.2 = c(4,5,6), only.pos = TRUE, max.cells.per.ident = 500, group.by = "orig.ident")

diffmark_0v1 <- FindMarkers(ARID1b_1M_rep1, ident.1 = 0, ident.2 = 2, only.pos = TRUE, max.cells.per.ident = 500, min.pct = 0.1)

#Graphical representations of cell types
FeaturePlot(ARID1b_1M_rep1, features = Thalamus[, "V1"], reduction = "tsne")

FeaturePlot(ARID1b_1M_rep1, features = oRG[, "V1"], reduction = "tsne")

FeaturePlot(ARID1b_1M_rep1, features = Cycling[, "V1"], reduction = "tsne")
cycling_genes <- c("MKI67", "TACC3", "TOP2A", "CDK1", "BIRC5", "CENPF")
FeaturePlot(ARID1b_1M_rep1, features = cycling_genes, reduction = "tsne")

FeaturePlot(ARID1b_1M_rep1, features = IPCs[, "V1"], reduction = "tsne")

FeaturePlot(ARID1b_1M_rep1, features = aRG[, "V1"], reduction = "tsne")
DotPlot(ARID1b_1M_rep1, features = aRG[, "V1"])
aRG_genes <- c("VIM", "HES1", "PAX6", "GLI3", "MKI67", "TACC3", "TOP2A", "ANXA1", "CDK1", "BIRC5", "CENPF")
FeaturePlot(ARID1b_1M_rep1, features = aRG_genes, reduction = "tsne")

differentiation_markers <- c("FOXG1", "SNAP25", "BCL11B", "SATB2", "GAD1", "DLX2", "VIM", "PAX6", "SPARCL1", "EMX1", "NKX2-1", "EOMES", "TOP2A")
FeaturePlot(ARID1b_1M_rep1, features = differentiation_markers, reduction = "tsne")

FeaturePlot(ARID1b_1M_rep1, features = Excitatory.Neurons[, "V1"], reduction = "tsne")
excitatory.neuron.genes <- c("NEUROD1", "NEUROD2", "NEUROD6", "NHLH1", "NHLH2", "TBR1", "GRIN2B", "GAP43", "ZEB2", "SNAP25")
FeaturePlot(ARID1b_1M_rep1, features = excitatory.neuron.genes, reduction = "tsne")

FeaturePlot(ARID1b_1M_rep1, features = INM[, "V1"], reduction = "tsne")     

DotPlot(ARID1b_1M_rep1, features = Cortical.Hem[, "V1"])
DoHeatmap(subset(ARID1b_1M_rep1, downsample = 100), features = Cortical.Hem[, "V1"]) + NoLegend()
FeaturePlot(ARID1b_1M_rep1, features = Cortical.Hem[, "V1"], reduction = "tsne")

FeaturePlot(ARID1b_1M_rep1, features = Excitatory.Neurons[,"V1"], reduction = "tsne")

FeaturePlot(ARID1b_1M_rep1, features = c("FEZF1, PCP4"), reduction = "tsne")

##RUN THIS ONE!! CONFIRM ARID1b PHENOTYPE
interneuron_genes <- c("DLX1", "DLX2", "GAD2")
FeaturePlot(ARID1b_1M_rep1, features = interneuron_genes, reduction = "tsne", split.by = "treat")
FeaturePlot(ARID1b_1M_rep1, features = IN[,"V1"], reduction = "tsne")

VlnPlot(ARID1b_1M_rep1, features = interneuron_genes, pt.size = 0, stack = T, sort=T, split.by = "treat")

differentiation_markers <- c("FOXG1", "SNAP25", "BCL11B", "SATB2", "GAD1", "DLX2", "VIM", "PAX6", "SPARCL1", "EMX1", "EOMES", "TOP2A", "NKX2-1")
FeaturePlot(ARID1b_1M_rep1, features = differentiation_markers, reduction = "tsne")
p1 <- VlnPlot(ARID1b_1M_rep1, features = differentiation_markers, pt.size = 0, stack = T)
p2 <- VlnPlot(Mito210_ARID1b_rep2, features = differentiation_markers, pt.size = 0, stack = T)

p1+p2

?VlnPlot
DotPlot(ARID1b_1M_rep1, features = "FOXG1")

cycling_genes <- c("MKI67", "TACC3", "TOP2A", "CDK1", "BIRC5", "CENPF")
FeaturePlot(ARID1b_1M_rep1, features = cycling_genes, reduction = "tsne")

cortical_hem_genes <- c("OTX2", "LMX1A", "WNT2B", "WNT5A", "WNT3A", "WLS", "FOXG1")
FeaturePlot(ARID1b_1M_rep1, features = cortical_hem_genes, reduction = "tsne")
?DotPlot()
DotPlot(ARID1b_1M_rep1, features = cortical_hem_genes)

IPC_genes <- c("EOMES", "PPP1R17", "PENK", "ELAVL4", "HES6", "NEUROD4")
FeaturePlot(ARID1b_1M_rep1, features = IPC_genes, reduction = "tsne")
DotPlot(ARID1b_1M_rep1, features = IPC_genes)

immature_neuronal_genes <- c("DCX", "NCAM1", "NEUROD1")
FeaturePlot(ARID1b_1M_rep1, features = immature_neuronal_genes, reduction = "tsne")

mature_neurons_genes <- c("ENO2", "RBFOX3", "MAP2", "TUBB3", "NEFL", "NEFM", "NEFH", "GAP43")
FeaturePlot(ARID1b_1M_rep1, features = mature_neurons_genes, reduction = "tsne")

synaptic_neuronal_genes <- c("DLG4", "SYP", "BSN")
FeaturePlot(ARID1b_1M_rep1, features = synaptic_neuronal_genes, reduction = "tsne")
FeaturePlot(ARID1b_1M_rep1, features = Synaptic[,"V1"], reduction = "tsne")
DotPlot(ARID1b_1M_rep1, features = synaptic_neuronal_genes)

FeaturePlot(ARID1b_1M_rep1, features = "CTIP2", reduction = "tsne")

FeaturePlot(ARID1b_1M_rep1, features = CajalRetzius[,"V1"], reduction = "tsne")

FeaturePlot(ARID1b_1M_rep1, features = RG[,"V1"], reduction = "tsne")

FeaturePlot(ARID1b_1M_rep1, features = CPN[,"V1"], reduction = "tsne")


####CLUSTER IDENTIFICATION SCRIPT####
#Data frame for relative prevalance of genes-related to specific cellular subtypes
genetic_marker_rarity <- data.frame()

#anti-hem
antihem.table <- data.frame()

for (i in 1:lengths(anti.hem)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["anti.hem", toString(i)] <- sum(antihem.table[,toString(i)]) / num.genes.antihem
}


# aRG
aRG.table <- data.frame()

for (i in 1:lengths(aRG)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["aRG", toString(i)] <- sum(aRG.table[,toString(i)]) / num.genes.aRG
}

#Astrocytes
astrocytes.table <- data.frame()

for (i in 1:lengths(Astrocytes)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Astrocytes", toString(i)] <- sum(astrocytes.table[,toString(i)]) / num.genes.astrocytes
}

#CFuPN_CPNs aspecific
CFuPN.CPN.table <- data.frame()

for (i in 1:lengths(CFUPN_CPN.aspecific)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["CFuPN.CPN", toString(i)] <- sum(CFuPN.CPN.table[,toString(i)]) / num.genes.CFuPN.CPN
}

#CFuPNs
CFuPN.table <- data.frame()

for (i in 1:lengths(CFuPNs)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["CFuPNs", toString(i)] <- sum(CFuPN.table[,toString(i)]) / num.genes.CFuPNs
}

#Cortical Hem
Cortical.hem.table <- data.frame()

for (i in 1:lengths(Cortical.Hem)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Cortical Hem", toString(i)] <- sum(CFuPN.CPN.table[,toString(i)]) / num.genes.cortical.hem
}


##Cycling progenitors
cycling.table <- data.frame()

for (i in 1:lengths(Cycling)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Cycling Progenitors", toString(i)] <- sum(cycling.table[,toString(i)]) / num.genes.cycling
}

#Ectoderm
ectoderm.table <- data.frame()

for (i in 1:lengths(Ectoderm)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Ectoderm", toString(i)] <- sum(ectoderm.table[,toString(i)]) / num.genes.ectoderm
}


##Excitatory neuronal progenitors
excitatory.neurons.table <- data.frame()

for (i in 1:lengths(Excitatory.Neurons)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Excitatory Neurons", toString(i)] <- sum(excitatory.neurons.table[,toString(i)]) / num.genes.excitatory.neurons
}

#Hypothalamus
hypothalamus.table <- data.frame()

for (i in 1:lengths(Hypothalamus)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Hypothalamus", toString(i)] <- sum(hypothalamus.table[,toString(i)]) / num.genes.hypothalamus
}

#INM
INM.table <- data.frame()

for (i in 1:lengths(INM)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["INM", toString(i)] <- sum(INM.table[,toString(i)]) / num.genes.INM
}

#IPC
IPC.table <- data.frame()

for (i in 1:lengths(IPCs)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["IPCs", toString(i)] <- sum(IPC.table[,toString(i)]) / num.genes.IPC
}

#Melanocytes
melanocyte.table <- data.frame()

for (i in 1:lengths(Melanocytes)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Melanocytes", toString(i)] <- sum(melanocyte.table[,toString(i)]) / num.genes.melanocytes
}

#Metabolic genes
metabolic.table <- data.frame()

for (i in 1:lengths(Metabolic.genes)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Metabolic Genes", toString(i)] <- sum(metabolic.table[,toString(i)]) / num.genes.metabolic
}

#Midbrain
midbrain.table <- data.frame()

for (i in 1:lengths(Midbrain)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Midbrain", toString(i)] <- sum(midbrain.table[,toString(i)]) / num.genes.midbrain
}

#Neural Crest
neural.crest.table <- data.frame()

for (i in 1:lengths(Neural.crest)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Neural Crest", toString(i)] <- sum(neural.crest.table[,toString(i)]) / num.genes.neural.crest
}

#Neuroepithelium
neuroepithelial.table <- data.frame()

for (i in 1:lengths(Neuroepithelial)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Neuroepithelial", toString(i)] <- sum(neuroepithelial.table[,toString(i)]) / num.genes.neuroeptihelial
}

#OLG
OLG.table <- data.frame()

for (i in 1:lengths(OLG)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["OLG", toString(i)] <- sum(OLG.table[,toString(i)]) / num.genes.OLG
}


#oRG
oRG.table <- data.frame()

for (i in 1:lengths(oRG)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["oRG", toString(i)] <- sum(oRG.table[,toString(i)]) / num.genes.oRG
}

#Retina
retina.table <- data.frame()

for (i in 1:lengths(Retina)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Retina", toString(i)] <- sum(retina.table[,toString(i)]) / num.genes.retina
}

#Striatum
striatum.table <- data.frame()

for (i in 1:lengths(Striatum)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Striatum", toString(i)] <- sum(striatum.table[,toString(i)]) / num.genes.striatum
}

#Subplate
subplate.table <- data.frame()

for (i in 1:lengths(Subplate)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Subplate", toString(i)] <- sum(subplate.table[,toString(i)]) / num.genes.subplate
}

#Synpatic
synaptic.table <- data.frame()

for (i in 1:lengths(Synaptic)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Synaptic", toString(i)] <- sum(synaptic.table[,toString(i)]) / num.genes.synaptic
}

#Thalamus
thalamus.table <- data.frame()

for (i in 1:lengths(Thalamus)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Thalamus", toString(i)] <- sum(thalamus.table[,toString(i)]) / num.genes.thalamus
}

#Ventral mix
ventral.mix.table <- data.frame()

for (i in 1:lengths(Ventral.mix)[1]) {
  for(k in 0:29){
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
for(i in 0:29){
  genetic_marker_rarity["Ventral Mix", toString(i)] <- sum(ventral.mix.table[,toString(i)]) / num.genes.ventral.mix
}

write.xlsx(genetic_marker_rarity, "~/Downloads/ARID1b_1M_rep1_3.15.21/ARID1b_1M_cluster_ident_frequency.xlsx", row.names = TRUE)










####Cluster assignment####

#To align Martina's cluster identities with my own. Correction was needed to have the cell names match
#the one's from Martina's Seurat object.
a <- substr(rownames(ARID1b_1M_rep1@meta.data), 1, nchar(rownames(ARID1b_1M_rep1@meta.data))-2)
ARID1b_1M_rep1$MartinaCellTypes = ARID1b_1M_rep1$CellType[a]

p1 = DimPlot(ARID1b_1M_rep1, reduction = "tsne", group.by = "MartinaCellTypes", label = TRUE)
p2 = DimPlot(ARID1b_1M_rep1, label = TRUE)

p1+p2

Idents(ARID1b_1M_rep1) <- "clusts_28PCs"
new.clusters.id <- c("Cycling aRG", "Cycling aRG", "Cycling aRG", "EN", "GABAergic", "EN", 
                     "Cycling aRG", "EN", "aRG", "EN", "Cortical Hem", "GABAergic", "Cortical Hem",
                     "EN", "Cycling GABAergic", "Cycling aRG", "GABAergic", "EN", "Cycling aRG", "GABAergic",
                     "Subcortical", "Subcortical EN", "Cycling aRG", "GABAergic", "Cycling IPC", "EN Cajal Retzius", 
                     "Cycling", "Cortical Hem", "aRG", "Cortical Hem")
names(new.clusters.id) <- levels(ARID1b_1M_rep1)
ARID1b_1M_rep1 <- RenameIdents(ARID1b_1M_rep1, new.clusters.id)

DimPlot(ARID1b_1M_rep1, reduction = "tsne", label = TRUE)
DimPlot(ARID1b_1M_rep1, reduction = "tsne", label = TRUE, split.by = "treat")

#cluster idents v2 --> Final cell typing 
v2.cluster.ids <- c("aRG", "aRG", "Cycling Progenitors", "Newborn PNs", "GABAergic Neurons 1",
                    "Newborn DL PNs", "Cycling Progenitors", "Newborn DL PNs", "aRG", "Newborn PNs",
                    "Cortical Hem", "GABAergic Neurons 1", "Cortical Hem", "IPC", "Cycling GABAergic Progenitors 1", "Subcortical", 
                    "GABAergic Progenitors 1", "Newborn DL PNs", "Cycling Progenitors", "GABAergic Neurons 1",
                    "Subcortical", "Subcortical", "Cycling Progenitors", "GABAergic Neurons 1",
                    "Cycling Progenitors", "Cajal-Retzius", "Cortical Hem", "Cortical Hem", "aRG", "Cortical Hem")
names(v2.cluster.ids) <- levels(ARID1b_1M_rep1)
ARID1b_1M_rep1 <- RenameIdents(ARID1b_1M_rep1, v2.cluster.ids)


ARID1b_1M_rep1$CellType = NA
#create a new column called "treatment" or whatever you want and set everything to NA at first
ARID1b_1M_rep1$CellType[ARID1b_1M_rep1$clusts_28PCs %in% c(0,1,8,28)] = "aRG"
ARID1b_1M_rep1$CellType[ARID1b_1M_rep1$clusts_28PCs %in% c(2,6,18,22,24)] = "Cycling Progenitors"
ARID1b_1M_rep1$CellType[ARID1b_1M_rep1$clusts_28PCs %in% c(3,9)] = "Newborn PNs"
ARID1b_1M_rep1$CellType[ARID1b_1M_rep1$clusts_28PCs %in% c(4,11,19,23)] = "GABAergic Neurons 1"
ARID1b_1M_rep1$CellType[ARID1b_1M_rep1$clusts_28PCs %in% c(5,7,17)] = "Newborn DL PNs"
ARID1b_1M_rep1$CellType[ARID1b_1M_rep1$clusts_28PCs %in% c(10,12,26,27,29)] = "Cortical Hem"
ARID1b_1M_rep1$CellType[ARID1b_1M_rep1$clusts_28PCs %in% c(13)] = "IPC"
ARID1b_1M_rep1$CellType[ARID1b_1M_rep1$clusts_28PCs %in% c(14)] = "Cycling GABAergic Progenitors 1"
ARID1b_1M_rep1$CellType[ARID1b_1M_rep1$clusts_28PCs %in% c(15,20,21)] = "Subcortical"
ARID1b_1M_rep1$CellType[ARID1b_1M_rep1$clusts_28PCs %in% c(16)] = "GABAergic Progenitors 1"
ARID1b_1M_rep1$CellType[ARID1b_1M_rep1$clusts_28PCs %in% c(25)] = "Cajal-Retzius"

DimPlot(ARID1b_1M_rep1, reduction = "tsne", group.by = "CellType", label = TRUE)


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

colsA1 = cols[match(levels(factor((ARID1b_1M_rep1)$CellType)),levels)]

DimPlot(ARID1b_1M_rep1, group.by = "CellType", cols=colsA1) + ggtitle('') + NoAxes() + theme_void()
DimPlot(ARID1b_1M_rep1, group.by = "CellType", cols=colsA1, split.by = "treat") + ggtitle("") + NoAxes()

saveRDS(ARID1b_1M_rep1, file = "clusteredSeur_ARID1b_1M.rds")

library(Seurat)
library(ggplot2)
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
ARID1b_1M_rep1 = readRDS(file = "clusteredSeur_ARID1B_1M.rds")

#Set Ident to cluster or celltypes, whatever groups you want to seperate before detecting DEGs in that group
Idents(ARID1b_1M_rep1) = "CellType"
id = "CellType"

#Set condition to the metadata column you want DEGs between, and base to the "wildtype" or base level of that column
condition = "treat"
base = "wt"

#set combineOn to the metadata column that contains the samples (i.e. different organoids)
combineOn = "orig.ident"

#This loop will run DE analysis for each cluster and save a .xlsx file for each!
for (id in levels(ARID1b_1M_rep1@active.ident)) {
  print(id)
  degs <- combineDE(ARID1b_1M_rep1, id=id, condition=condition, base=base, combineOn=combineOn)
  degs$gene = rownames(degs)
  if (length(degs)>0) {
    write_xlsx(degs, path=paste0(id,".DEGs.xlsx"), format_headers = T)
  }
}

##All GABAergic clusters were not considered because of their variable naming. 
##--> we may want to create an object with a combined GABAergic group

Idents(ARID1b_1M_rep1) = "clusts_28PCs"
id = "clusts_28PCs"
condition = "treat"
base = "wt"
combineOn = "orig.ident"

for (id in levels(ARID1b_1M_rep1@active.ident)) {
  print(id)
  degs <- combineDE(ARID1b_1M_rep1, id=id, condition=condition, base=base, combineOn=combineOn)
  degs$gene = rownames(degs)
  if (length(degs)>0) {
    write_xlsx(degs, path=paste0("clus",id,".DEGs.xlsx"), format_headers = T)
  }
}

##Since there was such a distinct wt/mt phenotype, many of the clusters did not get a deg list. 
##For example all of the GABAergic clusters had no noticeable cells in the WT orgs, so there was no deg generated to compare wt and mut groups


####GO and KEGG Pathway analysis####
library(dplyr)
library(Seurat)
library(ggplot2)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

#load all the libraries
library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(org.Hs.eg.db)
library(readr)
library(clusterProfiler)

#FindMarkers is ok, but it looks at cells as a single replicate --> use DESeq2

#Load in xlx file for CellType or Cluster of choice, subset by p-adjusted value <= 0.05
#make sure to load in the appropriate ARID1b 1M rep1!! DESeq2 generated files!!!!!!!!
#Isolate the genes names for GO/KEGG analyses

setwd("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_1M_scRNAseq_rep1_3.15.21")

#Read in DEGs for each and find overlaps and GO terms
allGenesUp = list()
allGenesDown = list()
datasets =  c("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_1M_scRNAseq_rep1_3.15.21/DESeq2 Cell Type",
              "~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_rep2_6.9.21/DESeq2 DEGs Cell Type")
names = c("ARID1b_1M_rep1", "ARID1b_1M_rep2")
dat = data.frame()
for (d in 1:length(datasets)) {
  dir=paste0(datasets[[d]],"DESeq2 Cell Type")
  xls = paste0(dir,".DEGs.xlsx")
  if (file.exists(xls)) {
    res = read_xlsx(xls)
    res = res[!is.na(res$pvalue),]
    res$dataset=names[[d]]
    print("ok")
    genesUp = res[res$padj<0.05 & res$log2FoldChange>0,"gene"]
    genesDown = res[res$padj<0.05 & res$log2FoldChange<0,"gene"]
    allGenesUp[[names[[d]]]] = genesUp$gene
    allGenesDown[[names[[d]]]] = genesDown$gene
    dat = rbind(dat,res)
  }
}


res = read_xlsx("Newborn DL PNs.DEGs.xlsx")
res = res[!is.na(res$pvalue),]
NewbornDLPN.genesUp = res[res$padj<0.05 & res$log2FoldChange>0,"gene"]
NewbornDLPN.genesUp <- unlist(NewbornDLPN.genesUp)
NewbornDLPN.genesDown = res[res$padj<0.05 & res$log2FoldChange<0,"gene"]
NewbornDLPN.genesDown <- unlist(NewbornDLPN.genesDown)

ego_NewbornDLPN_genesUp = enrichGO(gene = NewbornDLPN.genesUp, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
simplify_ego_NewbornDLPN_genesUp <- simplify(ego_NewbornDLPN_genesUp, cutoff = 0.7, by="p.adjust", select_fun=min)
dotplot(simplify_ego_NewbornDLPN_genesUp) + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle('GO ARID1b 1M rep1 Newborn DL PNs genesUp DEG')

ego_NewbornDLPN_genesDown = enrichGO(gene = NewbornDLPN.genesDown, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
simplify_ego_NewbornDLPN_genesDown <- simplify(ego_NewbornDLPN_genesDown, cutoff = 0.7, by="p.adjust", select_fun=min)
dotplot(simplify_ego_NewbornDLPN_genesDown) + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle('GO ARID1b 1M rep1 Newborn DL PNs genesDown DEG')

write.csv(data.frame(simplify_ego_NewbornDLPN_genesUp), file = "NewbornDLPN_genesUp_DEGs_GOterms.csv", row.names = FALSE)
write.csv(data.frame(simplify_ego_NewbornDLPN_genesDown), file = "NewbornDLPN_genesDown_DEGs_GOterms.csv", row.names = FALSE)

?dotplot

##aRG##
ARID1b_rep1_aRG_DEG_markers = subset(aRG_DEGs, padj <= 0.05)
ARID1b_rep1_aRG_DEG <- ARID1b_rep1_aRG_DEG_markers[,"gene"]
ARID1b_rep1_aRG_DEG <- unlist(ARID1b_rep1_aRG_DEG)

#subset now by positive and negative log fold change
aRG_DEG_posFC = subset(ARID1b_rep1_aRG_DEG_markers, log2FoldChange >= 0)
aRG_DEG_posFC_gene <- aRG_DEG_posFC[,"gene"]
aRG_DEG_posFC_gene <- unlist(aRG_DEG_posFC_gene)

aRG_DEG_negFC = subset(ARID1b_rep1_aRG_DEG_markers, log2FoldChange < 0)
aRG_DEG_negFC_gene <- aRG_DEG_negFC[,"gene"]
aRG_DEG_negFC_gene <- unlist(aRG_DEG_negFC_gene)

#aRG DEG general
ego_aRG = enrichGO(gene = ARID1b_rep1_aRG_DEG, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
simplify_ego_aRG <- simplify(ego_aRG, cutoff = 0.7, by="p.adjust", select_fun=min)
dotplot(ego_aRG) + theme(axis.text.x = element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle('GO ARID1b 1M rep1 aRG DEG')
dotplot(simplify_ego_aRG) + theme(axis.text.x = element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle('Simplified GO ARID1b 1M rep1 aRG DEG')
simplify_aRG_GOterm_list <- simplify_ego_aRG@result$Description

aRG_simple_summary <- data.frame(simplify_ego_aRG)
write.csv(aRG_simple_summary, file = "aRG_DEGs_GOterms.csv", row.names = FALSE)

#aRG DEG posFC mutatnt relative to wildtype
ego_aRG_posFC = enrichGO(gene = aRG_DEG_posFC_gene, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
dotplot(ego_aRG_posFC) + theme(axis.text.x = element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle('GO ARID1b 1M rep1 aRG DEG posFC')
simplify_ego_aRG_posFC <- simplify(ego_aRG_posFC, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(simplify_ego_aRG_posFC) + theme(axis.text.x = element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle('Simplified GO ARID1b 1M rep1 aRG DEG posFC')

aRG_posFC_summary <- data.frame(simplify_ego_aRG_posFC)
write.csv(aRG_posFC_summary, file = "aRG_posFC_DEGs_GOterms.csv", row.names = FALSE)

#aRG DEG negFC
ego_aRG_negFC = enrichGO(gene = aRG_DEG_negFC_gene, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
dotplot(ego_aRG_negFC) + theme(axis.text.x = element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle('GO ARID1b 1M rep1 aRG DEG negFC')
simplify_ego_aRG_negFC <- simplify(ego_aRG_negFC, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(simplify_ego_aRG_negFC) + theme(axis.text.x = element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle('Simplified GO ARID1b 1M rep1 aRG DEG negFC')
simplify_aRG_negFC_GOterm_list <- simplify_ego_aRG_negFC@result$Description

aRG_negFC_summary <- data.frame(simplify_ego_aRG_negFC)
write.csv(aRG_negFC_summary, file = "aRG_negFC_DEGs_GOterms.csv", row.names = FALSE)

rm(aRG_DEG_negFC, aRG_DEG_posFC, aRG_DEGs, ego_aRG, ego_aRG_negFC, ego_aRG_posFC, ARID1b_rep1_aRG_DEG_markers, simplify_aRG_negFC, simplify_aRG_posFC, simplify_aRG_posFC, simplify_ego_aRG, simplify_ego_aRG_negFC, simplify_ego_aRG_posFC, aRG_negFC_summary, aRG_simple_summary, aRG_posFC_summary)



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
dotplot(kk) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle('correlating genes\n no filter KEGG')

markers_mapIds_cluster4 <- mapIds(org.Hs.eg.db, cluster4_ARID1B_markers, 'ENTREZID', 'SYMBOL')
kk_cluster4 <- enrichKEGG(gene = markers_mapIds_cluster4, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_cluster4) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG aRG MT')


markers_mapIds_cluster9 <- mapIds(org.Hs.eg.db, cluster9_ARID1B_markers, 'ENTREZID', 'SYMBOL')
kk_cluster9 <- enrichKEGG(gene = markers_mapIds_cluster9, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_cluster4) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG aRG WT')

mapIds_WT_neuron <- mapIds(org.Hs.eg.db, WT_neuron, 'ENTREZID', 'SYMBOL')
kk_WT_neuron <- enrichKEGG(gene = mapIds_WT_neuron, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_WT_neuron) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG neuron WT')


mapIds_MT_neuron <- mapIds(org.Hs.eg.db, MT_neuron, 'ENTREZID', 'SYMBOL')
kk_MT_neuron <- enrichKEGG(gene = mapIds_MT_neuron, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_MT_neuron) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG neuron MT')

mapIds_neuron <- mapIds(org.Hs.eg.db, ARID1B_neuron, 'ENTREZID', 'SYMBOL')
kk_neuron <- enrichKEGG(gene = mapIds_neuron, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_neuron) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG ARID1B neuron')


mapIds_aRG <- mapIds(org.Hs.eg.db, ARID1B_aRG, 'ENTREZID', 'SYMBOL')
kk_aRG <- enrichKEGG(gene = mapIds_aRG, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_aRG) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG ARID1B aRG')


mapIds_diff_aRGs <- mapIds(org.Hs.eg.db, diffmark_aRG, 'ENTREZID', 'SYMBOL')
kk_diff_aRG <- enrichKEGG(gene = mapIds_diff_aRGs, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_diff_aRG) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG diff. markers aRG')

mapIds_diff_neuron <- mapIds(org.Hs.eg.db, diffmark_neuronal, 'ENTREZID', 'SYMBOL')
kk_diff_neuron <- enrichKEGG(gene = mapIds_diff_neuron, organism = "hsa", pvalueCutoff = 0.5)
dotplot(kk_diff_neuron) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle('correlating genes\n no filter KEGG diff. markers neuron')

#symbols <- as.matrix(symbols)
#symbols <- as.vector(symbols)


diffmark_aRG <- as.vector(diffmark_aRG)
diffmark_neuronal <- as.vector(diffmark_neuronal)




####Monocle3/Pseudotime Constructs####
#Load Monocle
install.packages("monocle3")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(Seurat)
library(cowplot)
library(ggplot2)
library(ggridges)

#Load data
setwd("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_1M_scRNAseq_rep1_3.15.21/")
ARID1b_1M_rep1 = readRDS("clusteredSeur_ARID1B_1M.rds")
fd = data.frame("gene_short_name" = rownames(ARID1b_1M_rep1))
rownames(fd) = rownames(ARID1b_1M_rep1)
cds <- new_cell_data_set(GetAssayData(ARID1b_1M_rep1,slot="counts"),
                         cell_metadata = ARID1b_1M_rep1@meta.data,
                         gene_metadata = fd)

#Subset so control cells equals mutant cells
table(colData(cds)$treat) #look at which one has fewer cells and downsample the other to that number
mutCells = rownames(colData(cds)[colData(cds)$treat=="mut",])
wtCells = sample(rownames(colData(cds)[colData(cds)$treat=="wt",]), 20552,replace = F)
cds = cds[,c(wtCells,mutCells)]


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
colData(cds)$clusts_28PCs = factor(colData(cds)$clusts_28PCs)
p2=plot_cells(cds, color_cells_by = "clusts_28PCs",show_trajectory_graph = F,
              label_groups_by_cluster = F,group_label_size = 4,cell_size = 0.5) + theme_void()
p3 = plot_cells(cds,show_trajectory_graph = T,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=1.5,
                cell_size = 0.5) + theme_void()
p4 = plot_cells(cds, color_cells_by = "treat", show_trajectory_graph = F, label_cell_groups = F,cell_size = 0.5) + theme_void()

p1
p2
p3
p4


?plot_cells()

colData(cds)$pseudotime <- pseudotime(cds)

df = as.data.frame(colData(cds))
df = df[is.finite(df$pseudotime),]

ggplot(df, aes(x=pseudotime, y=treat)) + 
  geom_density_ridges(aes(fill=treat)) +
  theme_void() +
  ggtitle('ARID1b 1M rep1 Pseudotime Original Split.by Treatment')

ggplot(df, aes(x=pseudotime, y=CellType)) + 
  geom_density_ridges(aes(fill=CellType)) +
  theme_void() + scale_color_manual(values=cols) +
  ggtitle('ARID1b 1M rep1 Pseudotime Original Split.by CellType')
##look strange with the altered pseudotime --> might be from a lack of connection to the main Cycling Prog/aRG body due to partition_qval=0.5 (from 1)

#subset out all non-cortical structures##
pseudo_subset_ARID1b_rep1 <- subset(ARID1b_1M_rep1, CellType %in% c("aRG", "Cajal-Retzius", "Cycling GABAergic Progenitors 1", "GABAergic Neurons 1", "GABAergic Progenitors 1", "Cycling Progenitors", "IPC", "Newborn PNs", "Newborn DL PNs"))
fd_subset_ARID1b_rep1 = data.frame("gene_short_name" = rownames(pseudo_subset_ARID1b_rep1))
rownames(fd_subset_ARID1b_rep1) = rownames(pseudo_subset_ARID1b_rep1)
cds_subset_ARID1b_rep1 <- new_cell_data_set(GetAssayData(pseudo_subset_ARID1b_rep1,slot="counts"),
                                cell_metadata = pseudo_subset_ARID1b_rep1@meta.data,
                                gene_metadata = fd_subset_ARID1b_rep1)

#Run Monocle
cds_subset_ARID1b_rep1 <- preprocess_cds(cds_subset_ARID1b_rep1, num_dim = 30, method = "PCA")

?preprocess_cds
#plot_pc_variance_explained(cds)
cds_subset_ARID1b_rep1 <- reduce_dimension(cds_subset_ARID1b_rep1,umap.fast_sgd=TRUE) #makes the UMAP
p1=plot_cells(cds_subset_ARID1b_rep1, color_cells_by = "CellType") #choose the metadata column you're interested in
cds_subset_ARID1b_rep1 <- cluster_cells(cds_subset_ARID1b_rep1, partition_qval=1) #playing with "partition_qval" from 0-1 will change how many different lineages your cells are split into
cds_subset_ARID1b_rep1 <- learn_graph(cds_subset_ARID1b_rep1)
cds_subset_ARID1b_rep1 = order_cells(cds_subset_ARID1b_rep1,root_pr_nodes=get_earliest_principal_node(cds_subset_ARID1b_rep1)) #Gives each cell a pseudotime value. Uses the helper function from above to set start point
p2=plot_cells(cds_subset_ARID1b_rep1,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=1.5)

p1
p2

##comparing WT and MT cells along the pseudotime trajectory via ridge plots
library(ggridges)
colData(cds_subset_ARID1b_rep1)$pseudotime <- pseudotime(cds_subset_ARID1b_rep1)
df = as.data.frame(colData(cds_subset_ARID1b_rep1))
df = df[is.finite(df$pseudotime),]
ggplot(df, aes(x=pseudotime, y=treat)) + 
  geom_density_ridges(aes(fill=treat, alpha=0.2)) +
  theme_classic() +
  ggtitle('ARID1b 1M rep1 Pseudotime Cortical Split.by Treatment')

ggplot(df, aes(x=pseudotime, y=CellType)) +
  geom_density_ridges(aes(fill=CellType)) +
  theme_classic() +
  ggtitle('ARID1b 1M rep1 Pseudotime Cortical Split.by CellType')


pseudo_MT_ARID1b_rep1 <- subset(ARID1b_1M_rep1, treat %in% c("mut"))
MT_ARID1b_rep1_subset <- subset(pseudo_MT_ARID1b_rep1, CellType %in% c("Newborn PNs", "Newborn DL PNs", "IPC", "GABAergic Progenitors 1", "GABAergic Neurons 1", "Cycling Progenitors", "Cycling GABAergic Progenitors 1", "Cajal-Retzius", "aRG"))

fd_MT_ARID1b_rep1 = data.frame("gene_short_name" = rownames(pseudo_MT_ARID1b_rep1))
rownames(fd_MT_ARID1b_rep1) = rownames(pseudo_MT_ARID1b_rep1)
cds_MT_ARID1b_rep1 <- new_cell_data_set(GetAssayData(pseudo_MT_ARID1b_rep1,slot="counts"),
                                            cell_metadata = pseudo_MT_ARID1b_rep1@meta.data,
                                            gene_metadata = fd_MT_ARID1b_rep1)

#Run Monocle
cds_MT_ARID1b_rep1 <- preprocess_cds(cds_MT_ARID1b_rep1, num_dim = 30, method = "PCA")

#plot_pc_variance_explained(cds)
cds_MT_ARID1b_rep1 <- reduce_dimension(cds_MT_ARID1b_rep1,umap.fast_sgd=TRUE) #makes the UMAP
cds_MT_ARID1b_rep1 <- cluster_cells(cds_MT_ARID1b_rep1, partition_qval=1) #playing with "partition_qval" from 0-1 will change how many different lineages your cells are split into
cds_MT_ARID1b_rep1 <- learn_graph(cds_MT_ARID1b_rep1)
p1=plot_cells(cds_MT_ARID1b_rep1, color_cells_by = "CellType") #choose the metadata column you're interested in
cds_MT_ARID1b_rep1 = order_cells(cds_MT_ARID1b_rep1,root_pr_nodes=get_earliest_principal_node(cds_MT_ARID1b_rep1)) #Gives each cell a pseudotime value. Uses the helper function from above to set start point
p2=plot_cells(cds_MT_ARID1b_rep1,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=1.5)

p1+p2

colData(cds_MT_ARID1b_rep1)$pseudotime <- pseudotime(cds_MT_ARID1b_rep1)
df = as.data.frame(colData(cds_MT_ARID1b_rep1))
df = df[is.finite(df$pseudotime),]

ggplot(df, aes(x=pseudotime, y=CellType)) +
  geom_density_ridges(aes(fill=CellType)) +
  theme_classic() +
  ggtitle('ARID1b 1M rep1 Pseudotime Cortical Split.by CellType')


fd_MT_ARID1b_rep1_subset = data.frame("gene_short_name" = rownames(MT_ARID1b_rep1_subset))
rownames(fd_MT_ARID1b_rep1_subset) = rownames(MT_ARID1b_rep1_subset)
cds_MT_ARID1b_rep1_subset <- new_cell_data_set(GetAssayData(MT_ARID1b_rep1_subset,slot="counts"),
                                        cell_metadata = MT_ARID1b_rep1_subset@meta.data,
                                        gene_metadata = fd_MT_ARID1b_rep1_subset)

#Run Monocle
cds_MT_ARID1b_rep1_subset <- preprocess_cds(cds_MT_ARID1b_rep1_subset, num_dim = 30, method = "PCA")

#plot_pc_variance_explained(cds)
cds_MT_ARID1b_rep1_subset <- reduce_dimension(cds_MT_ARID1b_rep1_subset,umap.fast_sgd=TRUE) #makes the UMAP
cds_MT_ARID1b_rep1_subset <- cluster_cells(cds_MT_ARID1b_rep1_subset, partition_qval=1) #playing with "partition_qval" from 0-1 will change how many different lineages your cells are split into
cds_MT_ARID1b_rep1_subset <- learn_graph(cds_MT_ARID1b_rep1_subset)
p1=plot_cells(cds_MT_ARID1b_rep1_subset, color_cells_by = "CellType") #choose the metadata column you're interested in
cds_MT_ARID1b_rep1_subset = order_cells(cds_MT_ARID1b_rep1_subset,root_pr_nodes=get_earliest_principal_node(cds_MT_ARID1b_rep1)) #Gives each cell a pseudotime value. Uses the helper function from above to set start point
p2=plot_cells(cds_MT_ARID1b_rep1_subset,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=1.5)

p1+p2

colData(cds_MT_ARID1b_rep1_subset)$pseudotime <- pseudotime(cds_MT_ARID1b_rep1_subset)
df = as.data.frame(colData(cds_MT_ARID1b_rep1_subset))
df = df[is.finite(df$pseudotime),]

ggplot(df, aes(x=pseudotime, y=CellType)) +
  geom_density_ridges(aes(fill=CellType)) +
  theme_classic() +
  ggtitle('ARID1b 1M rep1 Pseudotime MT Subset Cortical Split.by CellType')


##WT subset##
pseudo_WT_ARID1b_rep1 <- subset(ARID1b_1M_rep1, treat %in% c("wt"))
WT_ARID1b_rep1_subset <- subset(pseudo_WT_ARID1b_rep1, CellType %in% c("Newborn PNs", "Newborn DL PNs", "IPC", "Cycling Progenitors", "aRG", "GABAergic Neurons 1", "Cajal-Retzius"))

fd_WT_ARID1b_rep1 = data.frame("gene_short_name" = rownames(pseudo_WT_ARID1b_rep1))
rownames(fd_WT_ARID1b_rep1) = rownames(pseudo_WT_ARID1b_rep1)
cds_WT_ARID1b_rep1 <- new_cell_data_set(GetAssayData(pseudo_WT_ARID1b_rep1,slot="counts"),
                                        cell_metadata = pseudo_WT_ARID1b_rep1@meta.data,
                                        gene_metadata = fd_WT_ARID1b_rep1)

#Run Monocle
cds_WT_ARID1b_rep1 <- preprocess_cds(cds_WT_ARID1b_rep1, num_dim = 30, method = "PCA")

#plot_pc_variance_explained(cds)
cds_WT_ARID1b_rep1 <- reduce_dimension(cds_WT_ARID1b_rep1,umap.fast_sgd=TRUE) #makes the UMAP
cds_WT_ARID1b_rep1 <- cluster_cells(cds_WT_ARID1b_rep1, partition_qval=1) #playing with "partition_qval" from 0-1 will change how many different lineages your cells are split into
cds_WT_ARID1b_rep1 <- learn_graph(cds_WT_ARID1b_rep1)
p1=plot_cells(cds_WT_ARID1b_rep1, color_cells_by = "CellType") #choose the metadata column you're interested in
cds_WT_ARID1b_rep1 = order_cells(cds_WT_ARID1b_rep1,root_pr_nodes=get_earliest_principal_node(cds_WT_ARID1b_rep1)) #Gives each cell a pseudotime value. Uses the helper function from above to set start point
p2=plot_cells(cds_WT_ARID1b_rep1,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=1.5)

p1 + p2

colData(cds_WT_ARID1b_rep1)$pseudotime <- pseudotime(cds_WT_ARID1b_rep1)
df = as.data.frame(colData(cds_WT_ARID1b_rep1))
df = df[is.finite(df$pseudotime),]

ggplot(df, aes(x=pseudotime, y=CellType)) +
  geom_density_ridges(aes(fill=CellType)) +
  theme_classic() +
  ggtitle('ARID1b 1M rep1 Pseudotime WT Split.by CellType')


##Subset WT + cell type pseudo##
fd_WT_ARID1b_rep1_subset = data.frame("gene_short_name" = rownames(WT_ARID1b_rep1_subset))
rownames(fd_WT_ARID1b_rep1_subset) = rownames(WT_ARID1b_rep1_subset)
cds_WT_ARID1b_rep1_subset <- new_cell_data_set(GetAssayData(WT_ARID1b_rep1_subset,slot="counts"),
                                        cell_metadata = WT_ARID1b_rep1_subset@meta.data,
                                        gene_metadata = fd_WT_ARID1b_rep1)

#Run Monocle
cds_WT_ARID1b_rep1 <- preprocess_cds(cds_WT_ARID1b_rep1, num_dim = 30, method = "PCA")

#plot_pc_variance_explained(cds)
cds_WT_ARID1b_rep1 <- reduce_dimension(cds_WT_ARID1b_rep1,umap.fast_sgd=TRUE) #makes the UMAP
cds_WT_ARID1b_rep1 <- cluster_cells(cds_WT_ARID1b_rep1, partition_qval=1) #playing with "partition_qval" from 0-1 will change how many different lineages your cells are split into
cds_WT_ARID1b_rep1 <- learn_graph(cds_WT_ARID1b_rep1)
p1=plot_cells(cds_WT_ARID1b_rep1, color_cells_by = "CellType") #choose the metadata column you're interested in
cds_WT_ARID1b_rep1 = order_cells(cds_WT_ARID1b_rep1,root_pr_nodes=get_earliest_principal_node(cds_WT_ARID1b_rep1)) #Gives each cell a pseudotime value. Uses the helper function from above to set start point
p2=plot_cells(cds_WT_ARID1b_rep1,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=1.5)

p1 + p2

colData(cds_WT_ARID1b_rep1)$pseudotime <- pseudotime(cds_WT_ARID1b_rep1)
df = as.data.frame(colData(cds_WT_ARID1b_rep1))
df = df[is.finite(df$pseudotime),]

ggplot(df, aes(x=pseudotime, y=CellType)) +
  geom_density_ridges(aes(fill=CellType)) +
  theme_classic() +
  ggtitle('ARID1b 1M rep1 Pseudotime WT Split.by CellType')


pseudo_proliferatingsubset_ARID1b_rep1 <- subset(ARID1b_1M_rep1, CellType %in% c("aRG", "Cycling GABAergic Progenitors 1", "GABAergic Progenitors 1", "Cycling Progenitors", "IPC"))



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
cell_info <- ARID1b_1M_rep1
head(cell_info)

# Run propeller testing for cell type proportion differences between the two 
# groups
?propeller
propeller.test.ARID1b.rep1 <- propeller(clusters = ARID1b_1M_rep1$CellType, sample = ARID1b_1M_rep1$orig.ident, 
          group = ARID1b_1M_rep1$treat)
write.xlsx(propeller.test.ARID1b.rep1, file = "ARID1b_1M_rep1_propeller.xlsx")

propeller.test.ARID1b.combined <- propeller(clusters = ARID1b.1M.combined.obj$CellType, sample = ARID1b.1M.combined.obj$orig.ident, group = ARID1b.1M.combined.obj$treat)
plotCellTypeProps(clusters=ARID1b.1M.combined.obj$CellType, sample=ARID1b.1M.combined.obj$orig.ident) + theme_classic()



# Plot cell type proportions
plotCellTypeProps(clusters=ARID1b_1M_rep1$CellType, sample=ARID1b_1M_rep1$orig.ident) + theme_classic()
plotCellTypeProps(clusters=ARID1b_1M_rep1$CellType, sample=ARID1b_1M_rep1$treat) + theme_classic()


