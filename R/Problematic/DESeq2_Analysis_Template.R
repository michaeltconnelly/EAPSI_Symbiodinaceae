# Load in count table
library(SummarizedExperiment)
countdata <- assay(theTable)
coldata <- colData(theTable)
countdata <- read.delim("~/Documents/Research/EAPSI_Symbiodinaceae/outputs/SymC1_fC/heat_SymC1.counts.summary", header=FALSE)
head(countdata)
countdata <- colData(countdata)
rownames(coldata) <- coldata$run
colnames(countdata) <- coldata$run
head(coldata[, c("label1", "label2", "label3")])

# Construct data object
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ label1 + label2)
ddsFullCountTable

# Look for repeat runs of same sample
as.data.frame(colData(ddsFullCountTable)[, c("label1", "label2", "label3", "label4")])

# Collapse repeat runs of same sample
ddsCollapsed <- collapseReplicates(ddsFull, groupby = ddsFull$sample, run = ddsFull$run)
head(as.data.frame(colData(ddsCollapsed)[, c("sample", "runsCollapsed")]), 12)

# Confirm that counts for new object are equal to summed up ones
original <- rowSums(counts(ddsFullCountTable)[, ddsFullCountTable$label1 == "SRS308873"]) # (first replicate run)
all(original == counts(ddsCollapsed)[, "SRS308873"]) # if true will return TRUE

# Running the pipeline
dds <- ddsCollapsed[, ddsCollapsed$label4 == "48h"]
dds$time <- droplevels( dds$time ) # CONFUSED
dds$treatment <- relevel( dds$treatment, "Control" ) # CONFUSED
dds <- DESeq(dds)

# Inspect the results table
res <- results( dds )
res
mcols(res, use.names=TRUE)
res <- results( dds, contrast = c("treatment", "DPN", "Control") )
res

# MA plot (good for 2 group comparison)
plotMA( res, ylim = c(-1, 1) )

# Dispersion estimates plots
plotDispEsts( dds, ylim = c(1e-6, 1e1) )

# Histogram of p values
hist( res$pvalue, breaks=20, col="grey" )

# Basemean plot
qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
bins <- cut( res$baseMean, qs )
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")
attr(res,"filterThreshold")
plot(attr(res,"filterNumRej"),type="b",
     xlab="quantiles of 'baseMean'",
     ylab="number of rejections")

# Adding gene names
res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res$ensembl,
                  mart = ensembl )
idx <- match( res$ensembl, genemap$ensembl_gene_id )
res$entrez <- genemap$entrezgene[ idx ]
res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(res,4)

# Exporting results
res[1:2,]
write.csv( as.data.frame(res), file="results.csv" )

# Rlog transform
rld <- rlog( dds )
head( assay(rld) )
par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )
sampleDists <- dist( t( assay(rld) ) )
sampleDists
# (Heatmap)
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$treatment,
                                     rld$patient, sep="-" )
colnames(sampleDistMatrix) <- NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

# PCA
ramp <- 1:3/3
cols <- c(rgb(ramp, 0, 0),
          rgb(0, ramp, 0),
          rgb(0, 0, ramp),
          rgb(ramp, 0, ramp))
print( plotPCA( rld, intgroup = c( "patient", "treatment"), col=cols ) )

# Gene clustering
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
