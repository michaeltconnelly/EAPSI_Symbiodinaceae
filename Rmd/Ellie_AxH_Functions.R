#AxH_functions.R
#author: "Mike Connelly"
#date: "05/07/2020"

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("ComplexHeatmap")
library(tidyverse)
library(DESeq2)
library(geneplotter)
library(ggplot2)
library(grid)
library(pairheatmap)
library(pheatmap)

# Denoting res ---------------------------------------------------------------------------------------------

mycounts <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/outputs/DESeq-results/AntiandHeatCountsNew.csv")
metadata <- read_csv("~/Desktop/EAPSI_Symbiodinaceae/outputs/DESeq-results/AntiandHeatMeta.csv")
View(mycounts)
View(metadata)
class(mycounts)
class(metadata)
metadata$SampleID
names(mycounts)==metadata$SampleID
all(names(mycounts)==metadata$SampleID)
dds <- DESeqDataSetFromMatrix(countData = mycounts, colData = metadata, design = ~ Treatment + Colony)
dds <- DESeq(dds)
res <- results(dds)
res <- tbl_df(res)
View(res)

# EAPSI AXH Transcriptome Analysis Functions ---------------------------------------------------------------------------------
### PCA plot with custom PC axes----------------------------------------------------------------------------------
plotPCA.custom <-  function(object, intgroup="Treatment", ntop=500, returnData=FALSE, pcs = c(1,2))
{
  stopifnot(length(pcs) == 2)    ### added this to check number of PCs ####
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }
  # assemble the data for the plot
  ########## Here we just use the pcs object passed by the end user ####
  d <- data.frame(PC1=pca$x[,pcs[1]], PC2=pca$x[,pcs[2]], group=group, intgroup.df, name=colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  # extract loadings
}

### PCA plot formatted with aesthetics ---------------------------------------------------------------------------
ggPCA <- function(vsd, samples, condcolors, ntop = 500,  pclab = c(1,2)) {
  
  PCAtmtdata <- plotPCA.custom(vsd, intgroup = c("Colony", "Treatment"), ntop = 500, returnData = TRUE,  pcs = c(pclab[1],pclab[2]))
  #set factor orders 
  PCAtmtdata$Colony <- factor(PCAtmtdata$Colony, levels = c("HW1", "HW2", "WT1", "WT2"), ordered = TRUE)
  PCAtmtdata$Treatment <- factor(PCAtmtdata$Treatment, levels = c("control", "Heat", "Antibiotics", "Antibiotics.Heat"), ordered = TRUE)
  
  PCAtmtpercentVar <- round(100* attr(PCAtmtdata, "percentVar"))
  
  PCAplot <-  PCAtmtdata %>% ggplot(aes(PC1,PC2)) +
    geom_point(size=4, aes(fill = Treatment, shape = Colony), color = "black", stroke = 0.5, show.legend = TRUE) +
    xlab(paste0( "PC", pclab[1], ": ", PCAtmtpercentVar[pclab[1]], "% variance")) + 
    ylab(paste0( "PC", pclab[2], ": ", PCAtmtpercentVar[pclab[2]], "% variance")) + 
    coord_fixed(1) + 
    scale_fill_manual(values=condcolors, name="Treatment") + 
    scale_shape_manual(values=colshapes, name="Colony") +
    theme(legend.position = "right") +
    guides(fill = guide_legend(override.aes = list(fill = condcolors, shape = 21, alpha = 1, stroke = 0.5))) +
    ggtitle("Principal Component Analysis")
  
  PCAplot
}

PCA

### Function to plot color bar -------------------------------------------------------------------------------------
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }	
}

### Volcano plot for differential gene expression -----------------------------------------------------------------
volcanoplot <- function(res) {
  
  ##Highlight genes that have a padj < 0.05
  res$threshold <- ifelse(res$padj < 0.05 & res$log2FoldChange > 0, "Upregulated", ifelse(res$padj < 0.05 & res$log2FoldChange < 0, "Downregulated", "NA"))
  res$log10padj <- -log10(res$padj)
  dat_genes <- data.frame(cbind(res$log2FoldChange, res$log10padj, res$threshold), stringsAsFactors = FALSE)
  colnames(dat_genes) <- c("log2FoldChange", "log10padj", "threshold")
  row.names(dat_genes) <- res$IDGeneInfo
  #dat_genes <- dat_genes[order(dat_genes$log2FoldChange, decreasing = TRUE),]
  dat_genes$log2FoldChange <- as.numeric(dat_genes$log2FoldChange)
  dat_genes$log10padj <- as.numeric(dat_genes$log10padj)
  dat_genes$threshold <- factor(dat_genes$threshold, levels = c("Upregulated", "Downregulated", "NA"), ordered = TRUE)
  #Create volcanoplot
  gVolcano <- dat_genes %>% 
    ggplot(aes(log2FoldChange, log10padj)) + 
    geom_point(aes(color = threshold), alpha=0.7, size=2) +
    # scale_color_manual(values = DEGcolors) +
    scale_x_continuous(limits = c(-6,6), breaks = seq(-10,10,2)) + 
    ylim(c(0, 10)) +
    xlab("log2 fold change") +
    ylab("-log10 p-value") + 
    #geom_text_repel(data = dat_genes_LPS_ctrl[1:15, ], aes(label = rownames(dat_genes_LPS_ctrl[1:15, ])), color = "black", size = 2.5, box.padding = unit(0.35, "lines")) +
    theme_bw() +
    theme(legend.position = "none", 
          plot.title = element_text(size = 12, hjust = 0, vjust = 1),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10))
  print(gVolcano)
}

volcanoplot(res) # DEGcolors not defined...what scale color manual to use?

### Function to plot GO_MWU clustering ---------------------------------------------------------------------------

# From https://github.com/z0on/GO_MWU/blob/master/GO_MWU.R...

################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.

# Edit these to match your data file names: 
input=read_csv("~/Desktop/GO.csv") # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations=read_csv("~/Desktop/trial.txt") # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("~/Desktop/gomwu.functions.R")


# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# Plotting results
quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results

### Function to plot KOG delta ranks heatmap ---------------------------------------------------------------------
KOGheatmap <- function(deltaranks, pvals, ...) {
  deltaranks <- as.matrix(deltaranks)
  paletteLength <- 100
  myColor <- colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(paletteLength)
  myBreaks <- c(seq(min(deltaranks), 0, length.out = ceiling(paletteLength/2) + 1), 
                seq(max(deltaranks)/paletteLength, max(deltaranks), 
                    length.out = floor(paletteLength/2)))
  pheatmap(mat = deltaranks, display_numbers = pvals,
           cluster_cols = FALSE, cluster_rows = FALSE,
           treeheight_row = 15, treeheight_col = 15,
           border_color = "white", scale = "none",
           color = myColor, breaks = myBreaks, ...)
}

# EAPSI AXH Microbiome Analysis Functions ---------------------------------------------------------------------------------
### PCoA plot with spiders aesthetics

pcoa_spider <- function(ps, ord, centroids, type = "sample", axes = c(1,2)) {
  
  pcoa <- plot_ordination(ps, ord, type, axes) +
    #treatment ellipses
    #stat_ellipse(aes(color = Treatment, fill = Treatment), geom = "polygon", type = "norm", alpha = 0.0) + 
    #sample-centroid spiders paths
    geom_segment(data = centroids, mapping = aes(x = `Axis.1`, y = `Axis.2`, xend = c1, yend = c2),
                 lwd = 0.25, col = "dark grey") +
    #treatment centroid points
    geom_point(data = centroids, size = 4, aes(x = c1, y = c2, color = Treatment), fill = "black", shape = 21, stroke = 2, show.legend = TRUE) +
    #sample points
    geom_point(size = 3, aes(fill = Treatment, shape = Colony), color = "black", stroke = 0.5, show.legend = FALSE) +
    scale_color_manual(values = condcolors_AXH, name = "Treatment") +
    scale_fill_manual(values = condcolors_AXH, name = "Treatment") +
    scale_shape_manual(values = colshapes, name = "Colony") + 
    labs(x = "PC1", y = "PC2") +
    guides(shape = guide_legend(override.aes = list(shape = colshapes, alpha = 1, stroke = 0.5))) +
    theme(legend.spacing.y = unit(0, "cm"))
  return(pcoa)
}

