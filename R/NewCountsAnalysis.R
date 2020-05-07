# Loading in data

library(tidyverse)
library(DESeq2)
library(geneplotter)
library(ggplot2)
library(grid)
library(pairheatmap)
library(pheatmap)

mycounts <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatCountsNew.csv")
metadata <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatMeta.csv")
View(mycounts)
View(metadata)
class(mycounts)
class(metadata)
t.test(mycounts)
metadata$SampleID
names(mycounts)==metadata$SampleID
all(names(mycounts)==metadata$SampleID)
dds <- DESeqDataSetFromMatrix(countData = mycounts, colData = metadata, design = ~ Treatment + Colony)
dds <- DESeq(dds)
res <- results(dds)
res <- tbl_df(res)
View(res)
write_csv(res, "~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatResNew.csv")
labeledres <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatLabeledResNew.csv")

# Denoting significant expression as adjusted p value <= 0.05

labeledres <- labeledres %>% mutate(sig=padj<=0.05)
View(labeledres)
write_csv(labeledres, "~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatResforSigNew.csv")

# Reading in only significant data

sigres <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatSigResNew.csv")
View(sigres)

# Matching up significant scaffolds to gene models

scaffoldanno <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/ScaffoldAnno.csv")
View(scaffoldanno)
sigGOs <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatSigGOs.csv")
View(sigGOs)

# Matching up significant gene models to counts

sigGOcounts <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/SigGOCounts.csv")
View(sigGOcounts)

# Top 10 Up and Down Regulated

# Up

sigGOtop10up <- sigGOs %>% mutate(upreg=log2FoldChange>0)
View(sigGOtop10up)
sortup <- sigGOtop10up[order(sigGOtop10up[,15], decreasing = TRUE),]
View(sortup)
upreg <- sortup[1:4312,]
View(upreg[,13])
sortpadjup <- upreg[order(upreg[,13], decreasing = FALSE),]
View(sortpadjup[1:10,]) # These are the top 10 upregualted DEGs by padj, but because 6 are uncharacterized, let's keep going...
View(sortpadjup[1:16,]) # These include the top 10 known proteins
View(sortpadjup[1:60,]) # These include the top 10 proteins where the OS is symbiodinium (I just searched "symbiodinium" and selected those known GNs)

# Down

sigGOtop10down <- sigGOs %>% mutate(downreg=log2FoldChange<0)
View(sigGOtop10down)
sortdown <- sigGOtop10down[order(sigGOtop10down[,15], decreasing = TRUE),]
View(sortdown)
downreg <- sortdown[1:3766,]
View(downreg[,13])
sortpadjdown <- downreg[order(downreg[,13], decreasing = FALSE),]
View(sortpadjdown[1:10,]) # These are the top 10 downregualted DEGs by padj, but because 3 are uncharacterized, let's keep going...
View(sortpadjdown[1:15,]) # These include the top 10 known proteins
View(sortpadjdown[1:116,]) # These include the top 10 proteins where the OS is symbiodinium (I just searched "symbiodinium" and selected those known GNs)

# Generating MA plot

write_csv(labeledres, "~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/forMAmodNew.csv")
MAdata <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/forMANew.csv")
View(MAdata)
plotMA(MAdata)
plotMA(MAdata, ylim=c(-6,6))

# Generating PCA plot

StableSamples <- dds
rld <- rlogTransformation(StableSamples)
PCAmtdata <- plotPCA(rld, intgroup = c("Treatment"), returnData = TRUE)
View(PCAmtdata)
PCAmtpercentVar <- round(100* attr(PCAmtdata, "percentVar"))
condcolors <- c("#00CCCC", "#FFCC33", "#CC0033", "#0000FF", "#FF6600", "#FF66FF")
colshapes <- c(18, 9, 16, 10)
PCA <- ggplot(PCAmtdata, aes(PCAmtdata$PC1, PCAmtdata$PC2, color=StableSamples$Treatment, shape=StableSamples$Colony)) + geom_point(size=4, show.legend = TRUE) + xlab(paste0( "PC1: ", PCAmtpercentVar[1], "% variance")) + ylab(paste0( "PC2: ", PCAmtpercentVar[2], "% variance")) + coord_fixed() + ggtitle("Principal Component Analysis")
PCA + scale_color_manual(values=condcolors) + scale_shape_manual(values=colshapes)

# Heatmap by Counts

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)
selectdf <- as.data.frame(colData(dds)[,c("Colony", "Treatment")])
selectntd <- normTransform(dds)
pheatmap(assay(selectntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=selectdf, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)

# Heatmap by FoldChange

library(RColorBrewer)
library(genefilter)
library(gplots)

# By Treatment...

# (All)

all <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/NonRepFoldCounts.csv")
View(all)
allfolds <- all[3:15]
foldmatrix <- as.matrix(allfolds[, -1])
heatmap(foldmatrix)
breaksList = seq(-2, 2, by = 0.1)
pheatmap(foldmatrix, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = FALSE, annotation_legend = TRUE)

# (Top 20)

top20 <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Top20_FoldCounts_bypadj.csv")
top20 <- top20[3:15]
top20foldmatrix <- as.matrix(top20[, -1])
rownames(top20foldmatrix) <- top20$GN
heatmap(top20foldmatrix)
breaksList = seq(-5, 5, by = 0.3)
pheatmap(top20foldmatrix, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = FALSE, annotation_legend = TRUE)

# By Colony...

# (All)

bycol <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/FoldCounts_with_padj_byColony.csv")
bycol <- bycol[3:15]
bycolmatrix <- as.matrix(bycol[, -1])
pheatmap(bycolmatrix, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = FALSE, annotation_legend = TRUE)

# (Top 20)

rownames(bycolmatrix) <- bycol$GN
pheatmap(bycolmatrix[1:20, 1:12], color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = FALSE, annotation_legend = TRUE)

# Uncluster rows so treatment and colony can be compared

# (All)

pheatmap(foldmatrix, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = FALSE, cluster_rows = FALSE, annotation_legend = TRUE)
pheatmap(bycolmatrix, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = FALSE, cluster_rows = FALSE, annotation_legend = TRUE)

# (Top 20)

pheatmap(bycolmatrix[1:20, 1:12], color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = FALSE, cluster_rows = FALSE, annotation_legend = TRUE)
pheatmap(top20foldmatrix, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = FALSE, cluster_rows = FALSE, annotation_legend = TRUE)

# GO --- Work on this!

geneID2GOtab <- read_csv(file = "~/Desktop/Controlled.csv") 
View(geneID2GOtab)
geneID2GO <- geneID2GOtab[1:2]
View(geneID2GO)
geneUniverse <- names(geneID2GO) 
resultFisher <- runTest(geneID2GO, algorithm = "classic", statistic = "fisher")
