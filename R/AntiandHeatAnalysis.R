# Sorting Data

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
labeledres <- labeledres %>% mutate(sig=padj<0.05)
View(labeledres)
write_csv(labeledres, "~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatResforSig.csv")
sigres <- read_csv("~/Documents/Research/EAPSI_Symbiodinaceae/Ellie_R_CSVs/AntiandHeatSigRes.csv")
View(sigres)

