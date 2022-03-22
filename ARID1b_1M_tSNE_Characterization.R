p1=DimPlot(ARID1b_1M_rep1, reduction = "tsne", cols=colsA1, split.by = "treat") + ggtitle("\n") +theme_classic()
p2=DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", cols=colsA2, split.by = "treat") + ggtitle("\n") + theme_classic()

p3=DimPlot(ARID1b_1M_rep1, reduction = "tsne", cols=colsA1) + ggtitle("\n") +theme_classic() + NoLegend()
p4=DimPlot(Mito210_ARID1b_rep2, reduction = "tsne", cols=colsA2) + ggtitle("\n") +theme_classic() + NoLegend()




p1
p2

p3+p1

p4+p2


?DimPlot


