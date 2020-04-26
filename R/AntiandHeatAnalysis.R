# Loading in data

library(tidyverse)
library(DESeq2)
library(geneplotter)

mycounts <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatCounts.csv")
metadata <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatMeta.csv")
View(mycounts)
View(metadata)
class(mycounts)
class(metadata)
t.test(mycounts[-1])
metadata$SampleID
names(mycounts)[-1]==metadata$SampleID
all(names(mycounts)[-1]==metadata$SampleID)
dds <- DESeqDataSetFromMatrix(countData = mycounts[-1], colData = metadata, design = ~ Treatment + Colony)
dds <- DESeq(dds)
res <- results(dds)
res <- tbl_df(res)
View(res)
write_csv(res, "~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatRes.csv")
labeledres <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatLabeledRes.csv")

# Denoting significant expression as adjusted p value <= 0.05

labeledres <- labeledres %>% mutate(sig=padj<=0.05)
View(labeledres)
write_csv(labeledres, "~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatResforSig.csv")

# Reading in only significant data

sigres <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatSigRes.csv")
View(sigres)

# Matching up significant scaffolds to gene models

scaffoldanno <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/ScaffoldAnno.csv")
View(scaffoldanno)
sigproteins <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatSigProteins.csv")
View(sigproteins)

# Matching up significant gene models to counts

sigproteincounts <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatSigProteinsCounts.csv")
View(sigproteincounts)

# Generating MA plot

write_csv(labeledres, "~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/forMAmod.csv")
MAdata <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/forMA.csv")
View(MAdata)
plotMA(MAdata)
plotMA(MAdata, ylim=c(-6,6))

# Generating PCA plot

library(ggplot2)
StableSamples <- dds
rld <- rlogTransformation(StableSamples)
PCAmtdata <- plotPCA(rld, intgroup = c("Treatment"), returnData = TRUE)
PCAmtpercentVar <- round(100* attr(PCAmtdata, "percentVar"))
condcolors <- c("#00CCCC", "#FFCC33", "#CC0033", "#0000FF", "#FF6600", "#FF66FF")
colshapes <- c(18, 9, 16, 10)
PCA <- ggplot(PCAmtdata, aes(PCAmtdata$PC1, PCAmtdata$PC2, color=StableSamples$Treatment, shape=StableSamples$Colony)) + geom_point(size=4, show.legend = TRUE) + xlab(paste0( "PC1: ", PCAmtpercentVar[1], "% variance")) + ylab(paste0( "PC2: ", PCAmtpercentVar[2], "% variance")) + coord_fixed() + ggtitle("Principal Component Analysis")
PCA + scale_color_manual(values=condcolors) + scale_shape_manual(values=colshapes)

# Generating heat maps

library(grid)
library(pairheatmap)
library(pheatmap)

### Control (All Genotypes) 

gncontrol <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/GN_Control_Consolidate.txt", header=T)
View(gncontrol)
rownames(gncontrol) <- gncontrol$Gene
gncontrolmatrix <- as.matrix(gncontrol[, -1])
head(gncontrolmatrix)
heatmap(gncontrolmatrix, legend = TRUE)

### Heat (All Genotypes)

gnheat <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/GN_Heat_Consolidate.txt", header=T)
View(gnheat)
rownames(gnheat) <- gnheat$Gene
gnheatmatrix <- as.matrix(gnheat[, -1])
head(gnheatmatrix)
heatmap(gnheatmatrix)

### Anti (All Genotypes)

gnanti <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/GN_Anti_Consolidate.txt", header=T)
View(gnanti)
rownames(gnanti) <- gnanti$Gene
gnantimatrix <- as.matrix(gnanti[, -1])
head(gnantimatrix)
heatmap(gnantimatrix)

### Anti x Heat (All Genotypes)

gnantiheat <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/GN_AntiHeat_Consolidate.txt", header=T)
View(gnantiheat)
rownames(gnantiheat) <- gnantiheat$Gene
gnantiheatmatrix <- as.matrix(gnantiheat[, -1])
head(gnantiheatmatrix)
heatmap(gnantiheatmatrix)

###### Filtered; E Value = 0

### Control (All Genotypes) 

# Upregulated

upcontrol <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Control_Up_eval0.txt", header=T)
upcontrol <- upcontrol[1:12]
rownames(upcontrol) <- upcontrol$Gene
upcontrolmatrix <- as.matrix(upcontrol[, -1])
head(upcontrolmatrix)
heatmap(upcontrolmatrix)

# Downregulated

downcontrol <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Control_Down_eval0.txt", header=T)
downcontrol <- downcontrol[1:12]
rownames(downcontrol) <- downcontrol$Gene
downcontrolmatrix <- as.matrix(downcontrol[, -1])
head(downcontrolmatrix)
heatmap(downcontrolmatrix)

### Heat (All Genotypes)

# Upregulated

upheat <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Heat_Up_eval0.txt", header=T)
upheat <- upheat[1:13]
rownames(upheat) <- upheat$Gene
upheatmatrix <- as.matrix(upheat[, -1])
head(upheatmatrix)
heatmap(upheatmatrix)

# Downregulated

downheat <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Heat_Down_eval0.txt", header=T)
View(downheat)
downheat <- downheat[1:13]
rownames(downheat) <- downheat$Gene
downheatmatrix <- as.matrix(downheat[, -1])
head(downheatmatrix)
heatmap(downheatmatrix)

### Anti (All Genotypes)

# Upregulated

upanti <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Anti_Up_eval0.txt", header=T)
upanti <- upanti[1:13]
rownames(upanti) <- upanti$Gene
upantimatrix <- as.matrix(upanti[, -1])
head(upantimatrix)
heatmap(upantimatrix)

# Downregulated

downanti <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Anti_Down_eval0.txt", header=T)
downanti <- downanti[1:13]
rownames(downanti) <- downanti$Gene
downantimatrix <- as.matrix(downanti[, -1])
head(downantimatrix)
heatmap(downantimatrix)

### Anti x Heat (All Genotypes)

# Upregulated

upantiheat <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiHeat_Up_eval0.txt", header=T)
upantiheat <- upantiheat[1:13]
rownames(upantiheat) <- upantiheat$Gene
upantiheatmatrix <- as.matrix(upantiheat[, -1])
head(upantiheatmatrix)
heatmap(upantiheatmatrix)

# Downregulated

downantiheat <- read.table("~/Desktop/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiHeat_Down_eval0.txt", header=T)
downantiheat <- downantiheat[1:13]
rownames(downantiheat) <- downantiheat$Gene
downantiheatmatrix <- as.matrix(downantiheat[, -1])
head(downantiheatmatrix)
heatmap(downantiheatmatrix)



