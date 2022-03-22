

####ARID1B rep1####

#Percentage of cells in specific clusters per organoid
# install.packages("remotes")
remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
library(reshape)
library(reshape2)

Idents(ARID1b_1M_rep1) <- "CellType"


counts = as.matrix(table(ARID1b_1M_rep1$CellType,ARID1b_1M_rep1$orig.ident))
percent = t(t(counts)/colSums(counts))
percent = melt(percent,varnames = c('Cluster','Org'))
assigns = max.col(table(ARID1b_1M_rep1$orig.ident, ARID1b_1M_rep1$treat))
assigns = colnames(table(ARID1b_1M_rep1$orig.ident, ARID1b_1M_rep1$treat))[assigns]
names(assigns) = rownames(table(ARID1b_1M_rep1$orig.ident, ARID1b_1M_rep1$treat))
percent$treat = NA
for (o in 1:length(assigns)) { percent$treat[percent$Org==names(assigns)[[o]]] = assigns[[o]]}
percent$treat = factor(percent$treat, levels = c("wt","mut","ko"))
#percent$Org = factor(percent$Org, levels=c(5,6,1,2,3))#c("wt_1","wt_2","wt_3","mut_1","mut_2","mut_3"))
#percent$Cluster = relevel(percent$Cluster, ref="Cycling IN Progenitors")
#percent$Cluster = factor(as.character(percent$Cluster), levels=as.character(0:100))

levels =  c("Newborn DL PNs","Cajal-Retzius","Cortical Hem","Cycling Progenitors","IPC",
            "aRG","Newborn PNs","Subcortical", "Unknown", "Choroid Plexus/Cortical Hem",
            "CFuPNs","CPNs","oRG","PNs","GABAergic Neurons 2",
            "oRG/Astroglia","Astroglia","GABAergic Progenitors 2", "Cycling GABAergic Progenitors 2",
            "GABAergic Progenitors 1","Cycling GABAergic Progenitors 1","GABAergic Neurons 1")
cols = c("#f768a1","#ee8866","#bebada","#bbcc33","#fdb462",
         "#41ae76","#fa9fb5","#77aadd","darkgray", "powderblue",
         "#cc6677","#882255","#225522","#aa4499","#332288",
         "#009988","#0077bb","#5B65AE","#B87ACF",
         "turquoise","aquamarine","turquoise4")
cols = cols[match(levels(factor(percent$Cluster)),levels)]


percent.DLPN <- percent[percent$Cluster=="Newborn DL PNs",]

pdf("percent.Cells.in.all.CellTypes.pdf", height=2, width=9)
ggplot(percent[percent$Cluster!="Unknown",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values=cols)+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="Newborn DL PNs",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="#f768a1")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="Newborn PNs",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="#fa9fb5")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="GABAergic Neurons 1",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="turquoise4")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=90, hjust=0.5, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="Cycling GABAergic Progenitors 1",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="aquamarine")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=90, hjust=0.5, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="GABAergic Progenitors 1",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="turquoise")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=90, hjust=0.5, size=0), legend.position="none")

dev.off()


####ARID1B rep2####
Idents(Mito210_ARID1b_rep2) <- "CellType"


counts = as.matrix(table(Mito210_ARID1b_rep2$CellType,Mito210_ARID1b_rep2$orig.ident))
percent = t(t(counts)/colSums(counts))
percent = melt(percent,varnames = c('Cluster','Org'))
assigns = max.col(table(Mito210_ARID1b_rep2$orig.ident, Mito210_ARID1b_rep2$treat))
assigns = colnames(table(Mito210_ARID1b_rep2$orig.ident, Mito210_ARID1b_rep2$treat))[assigns]
names(assigns) = rownames(table(Mito210_ARID1b_rep2$orig.ident, Mito210_ARID1b_rep2$treat))
percent$treat = NA
for (o in 1:length(assigns)) { percent$treat[percent$Org==names(assigns)[[o]]] = assigns[[o]]}
percent$treat = factor(percent$treat, levels = c("wt","mut","ko"))
#percent$Org = factor(percent$Org, levels=c(5,6,1,2,3))#c("wt_1","wt_2","wt_3","mut_1","mut_2","mut_3"))
#percent$Cluster = relevel(percent$Cluster, ref="Cycling IN Progenitors")
#percent$Cluster = factor(as.character(percent$Cluster), levels=as.character(0:100))

levels =  c("Newborn DL PNs","Cajal-Retzius","Cortical Hem","Cycling Progenitors","IPC",
            "aRG","Newborn PNs","Subcortical", "Unknown", "Choroid Plexus/Cortical Hem",
            "CFuPNs","CPNs","oRG","PNs","GABAergic Neurons 2",
            "oRG/Astroglia","Astroglia","GABAergic Progenitors 2", "Cycling GABAergic Progenitors 2",
            "GABAergic Progenitors 1","Cycling GABAergic Progenitors 1","GABAergic Neurons 1")
cols = c("#f768a1","#ee8866","#bebada","#bbcc33","#fdb462",
         "#41ae76","#fa9fb5","#77aadd","darkgray", "powderblue",
         "#cc6677","#882255","#225522","#aa4499","#332288",
         "#009988","#0077bb","#5B65AE","#B87ACF",
         "turquoise","aquamarine","turquoise4")
cols = cols[match(levels(factor(percent$Cluster)),levels)]


percent.DLPN <- percent[percent$Cluster=="Newborn DL PNs",]

ggplot(percent[percent$Cluster!="Unknown",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values=cols)+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="Newborn DL PNs",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="#f768a1")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="Newborn PNs",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="#fa9fb5")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="GABAergic Neurons 1",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="turquoise4")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=90, hjust=0.5, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="Cycling GABAergic Progenitors 1",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="aquamarine")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=90, hjust=0.5, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="GABAergic Progenitors 1",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="turquoise")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=90, hjust=0.5, size=0), legend.position="none")

dev.off()


####PTEN rep1####

counts = as.matrix(table(PTEN.rep1.obj$CellType,PTEN.rep1.obj$orig.ident))
percent = t(t(counts)/colSums(counts))
percent = melt(percent,varnames = c('Cluster','Org'))
assigns = max.col(table(PTEN.rep1.obj$orig.ident, PTEN.rep1.obj$treat))
assigns = colnames(table(PTEN.rep1.obj$orig.ident, PTEN.rep1.obj$treat))[assigns]
names(assigns) = rownames(table(PTEN.rep1.obj$orig.ident, PTEN.rep1.obj$treat))
percent$treat = NA
for (o in 1:length(assigns)) { percent$treat[percent$Org==names(assigns)[[o]]] = assigns[[o]]}
percent$treat = factor(percent$treat, levels = c("wt","mt","ko"))
#reorder the orgs in wt --> mutant order
percent$Org = factor(percent$Org, levels=c(4,5,1,2,3))#c("wt_1","wt_2","wt_3","mut_1","mut_2","mut_3"))
#percent$Cluster = relevel(percent$Cluster, ref="Cycling IN Progenitors")
#percent$Cluster = factor(as.character(percent$Cluster), levels=as.character(0:100))

levels =  c("Newborn DL PNs","Cajal-Retzius","Cortical Hem","Cycling Progenitors","IPC",
            "aRG","Newborn PNs","Subcortical", "Unknown", "Choroid Plexus/Cortical Hem",
            "CFuPNs","CPNs","oRG","PNs","GABAergic Neurons 2",
            "oRG/Astroglia","Astroglia","GABAergic Progenitors 2", "Cycling GABAergic Progenitors 2",
            "GABAergic Progenitors 1","Cycling GABAergic Progenitors 1","GABAergic Neurons 1")
cols = c("#f768a1","#ee8866","#bebada","#bbcc33","#fdb462",
         "#41ae76","#fa9fb5","#77aadd","darkgray", "powderblue",
         "#cc6677","#882255","#225522","#aa4499","#332288",
         "#009988","#0077bb","#5B65AE","#B87ACF",
         "turquoise","aquamarine","turquoise4")
cols = cols[match(levels(factor(percent$Cluster)),levels)]


percent.DLPN <- percent[percent$Cluster=="Newborn DL PNs",]

pdf("percent.Cells.in.all.CellTypes.pdf", height=2, width=9)
ggplot(percent[percent$Cluster!="Unknown",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values=cols)+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=14), legend.position="none")

ggplot(percent[percent$Cluster=="Newborn DL PNs",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="#f768a1")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="Newborn PNs",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="#fa9fb5")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=0.5, size=14), legend.position="none")

dev.off()

####PTEN rep2####
counts = as.matrix(table(PTEN.rep2.obj$CellType,PTEN.rep2.obj$orig.ident))
percent = t(t(counts)/colSums(counts))
percent = melt(percent,varnames = c('Cluster','Org'))
assigns = max.col(table(PTEN.rep2.obj$orig.ident, PTEN.rep2.obj$treat))
assigns = colnames(table(PTEN.rep2.obj$orig.ident, PTEN.rep2.obj$treat))[assigns]
names(assigns) = rownames(table(PTEN.rep2.obj$orig.ident, PTEN.rep2.obj$treat))
percent$treat = NA
for (o in 1:length(assigns)) { percent$treat[percent$Org==names(assigns)[[o]]] = assigns[[o]]}
percent$treat = factor(percent$treat, levels = c("wt","mt","ko"))
#reorder the orgs in wt --> mutant order
percent$Org = factor(percent$Org, levels=c(4,5,1,2,3))#c("wt_1","wt_2","wt_3","mut_1","mut_2","mut_3"))
#percent$Cluster = relevel(percent$Cluster, ref="Cycling IN Progenitors")
#percent$Cluster = factor(as.character(percent$Cluster), levels=as.character(0:100))

levels =  c("Newborn DL PNs","Cajal-Retzius","Cortical Hem","Cycling Progenitors","IPC",
            "aRG","Newborn PNs","Subcortical", "Unknown", "Choroid Plexus/Cortical Hem",
            "CFuPNs","CPNs","oRG","PNs","GABAergic Neurons 2",
            "oRG/Astroglia","Astroglia","GABAergic Progenitors 2", "Cycling GABAergic Progenitors 2",
            "GABAergic Progenitors 1","Cycling GABAergic Progenitors 1","GABAergic Neurons 1")
cols = c("#f768a1","#ee8866","#bebada","#bbcc33","#fdb462",
         "#41ae76","#fa9fb5","#77aadd","darkgray", "powderblue",
         "#cc6677","#882255","#225522","#aa4499","#332288",
         "#009988","#0077bb","#5B65AE","#B87ACF",
         "turquoise","aquamarine","turquoise4")
cols = cols[match(levels(factor(percent$Cluster)),levels)]


percent.DLPN <- percent[percent$Cluster=="Newborn DL PNs",]

pdf("percent.Cells.in.all.CellTypes.pdf", height=2, width=9)
ggplot(percent[percent$Cluster!="Unknown",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values=cols)+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=14), legend.position="none")

ggplot(percent[percent$Cluster=="Newborn DL PNs",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="#f768a1")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=0), legend.position="none")

ggplot(percent[percent$Cluster=="Newborn PNs",],aes(y=value, x=Org,fill=Cluster,color=treat)) + 
  geom_bar(stat="identity",show.legend=T,size=0.8)+
  scale_fill_manual(values="#fa9fb5")+
  #geom_col_pattern(aes(pattern_alpha=treat),pattern="crosshatch",pattern_color="NA",color="NA",pattern_type="tile")+
  #scale_pattern_alpha_manual(values=c(0,1))+
  scale_color_manual(values=c("darkgray","black"))+
  facet_grid(. ~ Cluster) + theme_classic() +
  labs(y="Percentage of Cells", x = "Organoid") +
  theme(axis.text.x = element_text(angle=70, hjust=0.5, size=14), legend.position="none")

dev.off()


