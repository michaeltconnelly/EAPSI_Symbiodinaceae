# Loading in data

library(tidyverse)
library(DESeq2)
library(geneplotter)

mycounts <- read_csv("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatCounts.csv")
metadata <- read_csv("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatMeta.csv")
View(mycounts)
View(metadata)
class(mycounts)
class(metadata)
metadata$SampleID
names(mycounts)[-1]==metadata$SampleID
all(names(mycounts)[-1]==metadata$SampleID)
dds <- DESeqDataSetFromMatrix(countData = mycounts[-1], colData = metadata, design = ~ Treatment + Colony)
dds <- DESeq(dds)
res <- results(dds)
res <- tbl_df(res)
View(res)
write_csv(res, "~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatRes.csv")
labeledres <- read_csv("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatLabeledRes.csv")

# Denoting significant expression as adjusted p value <= 0.05

labeledres <- labeledres %>% mutate(sig=padj<=0.05)
View(labeledres)
write_csv(labeledres, "~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatResforSig.csv")

# Reading in only significant data

sigres <- read_csv("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatSigRes.csv")
View(sigres)

# Matching up significant scaffolds to proteins

scaffoldanno <- read_csv("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/ScaffoldAnno.csv")
View(scaffoldanno)
sigproteins <- read_csv("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatSigProteins.csv")
View(sigproteins)

# Matching up significant proteins to counts

sigproteincounts <- read_csv("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatSigProteinsCounts.csv")
View(sigproteincounts)

# Matching up significant protein counts to treatments - TBC

metadata$SampleID
ls(sigproteincounts[15:61])
all(ls(sigproteincounts[15:61])==metadata$SampleID)

# GENERATING HEAT MAP (All Unique DEGs)

library(grid)
library(pairheatmap)
library(pheatmap)
help(pheatmap)

### Control (All Genotypes)

gncontrol <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/GN_Control_Consolidate.txt", header=T)
View(gncontrol)
rownames(gncontrol) <- gncontrol$Gene
gncontrolmatrix <- as.matrix(gncontrol[, -1])
head(gncontrolmatrix)
heatmap(gncontrolmatrix, legend = TRUE)

### Heat (All Genotypes)

gnheat <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/GN_Heat_Consolidate.txt", header=T)
View(gnheat)
rownames(gnheat) <- gnheat$Gene
gnheatmatrix <- as.matrix(gnheat[, -1])
head(gnheatmatrix)
heatmap(gnheatmatrix)

### Anti (All Genotypes)

gnanti <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/GN_Anti_Consolidate.txt", header=T)
View(gnanti)
rownames(gnanti) <- gnanti$Gene
gnantimatrix <- as.matrix(gnanti[, -1])
head(gnantimatrix)
heatmap(gnantimatrix)

### Anti x Heat (All Genotypes)

gnantiheat <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/GN_AntiHeat_Consolidate.txt", header=T)
View(gnantiheat)
rownames(gnantiheat) <- gnantiheat$Gene
gnantiheatmatrix <- as.matrix(gnantiheat[, -1])
head(gnantiheatmatrix)
heatmap(gnantiheatmatrix)

###### Filtered; E Value = 0

### Control (All Genotypes)

# Upregulated

upcontrol <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Control_Up_eval0.txt", header=T)
upcontrol <- upcontrol[1:12]
rownames(upcontrol) <- upcontrol$Gene
upcontrolmatrix <- as.matrix(upcontrol[, -1])
head(upcontrolmatrix)
heatmap(upcontrolmatrix)

# Downregulated

downcontrol <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Control_Down_eval0.txt", header=T)
downcontrol <- downcontrol[1:12]
rownames(downcontrol) <- downcontrol$Gene
downcontrolmatrix <- as.matrix(downcontrol[, -1])
head(downcontrolmatrix)
heatmap(downcontrolmatrix)

### Heat (All Genotypes)

# Upregulated

upheat <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Heat_Up_eval0.txt", header=T)
upheat <- upheat[1:13]
rownames(upheat) <- upheat$Gene
upheatmatrix <- as.matrix(upheat[, -1])
head(upheatmatrix)
heatmap(upheatmatrix)

# Downregulated

downheat <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Heat_Down_eval0.txt", header=T)
View(downheat)
downheat <- downheat[1:13]
rownames(downheat) <- downheat$Gene
downheatmatrix <- as.matrix(downheat[, -1])
head(downheatmatrix)
heatmap(downheatmatrix)

### Anti (All Genotypes)

# Upregulated

upanti <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Anti_Up_eval0.txt", header=T)
upanti <- upanti[1:13]
rownames(upanti) <- upanti$Gene
upantimatrix <- as.matrix(upanti[, -1])
head(upantimatrix)
heatmap(upantimatrix)

# Downregulated

downanti <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Anti_Down_eval0.txt", header=T)
downanti <- downanti[1:13]
rownames(downanti) <- downanti$Gene
downantimatrix <- as.matrix(downanti[, -1])
head(downantimatrix)
heatmap(downantimatrix)

### Anti x Heat (All Genotypes)

# Upregulated

upantiheat <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiHeat_Up_eval0.txt", header=T)
upantiheat <- upantiheat[1:13]
rownames(upantiheat) <- upantiheat$Gene
upantiheatmatrix <- as.matrix(upantiheat[, -1])
head(upantiheatmatrix)
heatmap(upantiheatmatrix)

# Downregulated

downantiheat <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiHeat_Down_eval0.txt", header=T)
downantiheat <- downantiheat[1:13]
rownames(downantiheat) <- downantiheat$Gene
downantiheatmatrix <- as.matrix(downantiheat[, -1])
head(downantiheatmatrix)
heatmap(downantiheatmatrix)

######### Filtered; E Value = 0

###### Anti x Heat

# Hw1

library(pheatmap)

antiheathw1 <- read.table("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/Hw1_AntiHeat_20most.txt", header=T)
View(antiheathw1[1:4])
# DEG <- antiheathw1[5:6]
# View(DEG)
antiheathw1 <- antiheathw1[1:4]
rownames(antiheathw1) <- antiheathw1$Gene
antiheathw1matrix <- as.matrix(antiheathw1[, -1])
pheatmap(antiheathw1matrix)
# anno_colors <- list(DEG=c(Upregulated="red", Downregulated="blue"))
# pheatmap(antiheathw1matrix, annotation_colors = anno_colors)
# pheatmap(antiheathw1matrix)

