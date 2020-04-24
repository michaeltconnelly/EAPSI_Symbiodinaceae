# Sorting Data


library(tidyverse)
library(DESeq2)
library(geneplotter)
mycounts <- read_csv("~/Documents/Research/alltreatments_SymC1_IDs.csv")
metadata <- read_csv("~/Documents/Research/EAPSIsamples_stable_meta.csv")
anno <- read_csv("~/Documents/Research/heat_SymC1.csv")
View(mycounts)
View(metadata)
class(mycounts)
class(metadata)
metadata$SampleID
names(mycounts)[-1]==metadata$SampleID
all(names(mycounts)[-1]==metadata$SampleID)
dds <- DESeqDataSetFromMatrix(countData = mycounts[-1], colData = metadata, design = ~ Treatment)
dds <- DESeq(dds)
res <- results(dds)
res <- tbl_df(res)
View(res)
write_csv(res, "sigres.csv")
geneids <- read_csv("~/Documents/Research/GeneIDs.csv")
View(geneids)
# merge(res, geneids) -- Memory exhausted; let's try creating a CSV
labeledres <- read_csv("~/Documents/Research/labeled_results.csv")
View(labeledres)
pvaluesortres <- labeledres %>% arrange(pvalue)
View(pvaluesortres)
pvaluesortres %>% filter(padj<=0.05) %>% write_csv("signif.csv") # Only 9 significant genes?!
labeledres <- labeledres %>% mutate(sig=padj<0.05)
labeledres %>% group_by(sig) %>% summarize(n=n())
View(labeledres)


# MA

'''
write_csv(labeledres, "forMAmod.csv")
# MAdata <- read_csv("~/Documents/Research/MAdata.csv")
MAdata <- read_csv("~/Documents/Research/tryok.csv")
View(MAdata)
# plotMA(MAdata, ylim=c(-1,1), xlim=c(50,1000)) # Modify format a bit now
plotMA(MAdata, ylim = c(-2,2), colNonSig = "gray", colSig = "red", colLine = "blue", log = "x", cex=0.45, xlab="mean expression", ylab="log fold change")
'''

library("tidyverse")
library("DESeq2")
library("tximport")
library("readxl")


View(metadata)


dds <- DESeqDataSetFromMatrix(mycounts[-1],
                                colData = metadata,
                                design = ~Batch + Colony + Treatment + Colony*Treatment)
dds$Treatment <- factor(dds$Treatment, levels = c("control", "Antibiotics", "Heat", "Antibiotics.Heat"))
relevel(dds$Treatment, ref = "control")
keep <- rowSums(counts(dds), na.rm = TRUE) >= 10
dds <- dds[keep,]
dds %>% na.exclude
RemovedNAs <- na.exclude(dds$Treatment)
View(RemovedNAs)


DESeqDataSetFromMatrix(RemovedNAs, colData=metadata)
RemovedNAs <- DESeq(RemovedNAs)
res <- results(dds)


# PCA

library(EDASeq)
library(ggplot2)

# PCAmat <- as.matrix(metadata)
# PCAdata <- newSeqExpressionSet(PCAmat,treatData=AnnotatedDataFrame(data.frame(Treatment=factor(c("Heat", "Antibiotics.Heat", "Antibiotics")), row.names=colnames(mycounts[-1]))))
# plotPCA(object, intgroup=c("condition", "type"))
rld <- rlogTransformation(dds)
# plotPCA(rld, intgroup = c("Treatment", "Colony"), returnData = FALSE) # Changed to FALSE
PCAtmtdata <- plotPCA(rld, intgroup = c("Treatment", "Colony"), returnData = TRUE)
View(PCAtmtdata)
# cols<-!(colnames(PCAtmtdata) %in% c("group","name"))
# CondensedData <-subset(PCAtmtdata,,cols)
# View(CondensedData)
PCAtmtpercentVar <- round(100* attr(PCAtmtdata, "percentVar"))
condcolors <- c("#00CCCC", "#FFCC33", "#CC0033", "#0000FF", "#FF6600", "#FF66FF")
colshapes <- c(18, 9, 16, 10)
PCA <- ggplot(PCAtmtdata, aes(PCAtmtdata$PC1, PCAtmtdata$PC2)) + geom_point(size=4, show.legend = TRUE) + xlab(paste0( "PC1: ", PCAtmtpercentVar[1], "% variance")) + ylab(paste0( "PC2: ", PCAtmtpercentVar[2], "% variance")) + coord_fixed() + ggtitle("Principal Component Analysis", subtitle = "Overall Symbiodinium transcriptome expression data")
PCA + scale_color_manual(values=condcolors) + scale_shape_manual(values=colshapes)
# habillage=iris$Treatment
# PCA + habillage
trycolors <- ggplot(PCAtmtdata, aes(x = PC1, y = PC2, color = factor(Treatment), shape = factor(Colony))) + 
  geom_point(size =3, aes(fill=factor(Colony))) +
  geom_point(size =3) + 
  scale_shape_manual(values=c(21,22)) + 
  scale_alpha_manual(values=c("F"=0, "M"=1))
trycolors
