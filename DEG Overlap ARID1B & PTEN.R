#PTEN and ARID1B 1M
#Read in DEGs for each and find overlaps and GO terms
install.packages("corrplot")
library(corrplot)

allGenes = list()
dat = data.frame()
folders =  c("~/Documents/Harvard/Arlotta Lab/Mito210_ARID1b_1M_Harmony/DESeq2 Cell Type/",
             "~/Documents/Harvard/Arlotta Lab/Mito210_PTEN_1M_Harmony/DESeq2 Cell Type/")
names = c("ARID1B 1M","PTEN 1M")

for(i in 1:length(names)) {
  files = list.files(path=folders[[i]], full.names = F)
  for (d in 1:length(files)) {
    ct = gsub(".DEGs.xlsx"," ",files[[d]])
    xls = paste0(folders[[i]],files[[d]])
    res = read_xlsx(xls)
    res$dataset=paste0(names[[i]],"-",ct)
    dat = rbind(dat,res)
    genes = res[!is.na(res$padj) & res$padj<0.05,"gene"]
    allGenes[[paste0(names[[i]],"-",ct)]] = genes$gene
  }
}


#Get genes that overlap
overlap = sapply(allGenes, function(x) sapply(allGenes, function(y) sum(y %in% x)))

#hypergeometric
ps = data.frame()
for (r in rownames(overlap)) {
  for (c in colnames(overlap)) {
    universe = length(intersect(dat[dat$dataset==r,]$gene, dat[dat$dataset==c,]$gene))
    rowsize = length(allGenes[[r]])
    colsize = length(allGenes[[c]])
    numoverlap = overlap[r,c]
    p = phyper(numoverlap-1, rowsize, universe-rowsize, colsize, lower.tail = F)
    ps[r,c] = p
  }
}
numTests = sum(1:(nrow(overlap)-1))
padjs = apply(ps, c(1,2), function (x) min(x*numTests, 1))
padjs = apply(padjs, c(1,2), function (x) max(x, 10^-100))
logps = as.matrix(-1* log10(padjs))

pdf("celltypes-SUVandARID1m-allDEGs-Overlap-maxp.pdf")
corrplot(logps, is.corr=F, type = "lower", diag = F, tl.col = "black", cl.length=2)
dev.off()