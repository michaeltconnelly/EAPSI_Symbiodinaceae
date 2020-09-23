# this is the final script in GO analysis using MWU test, as in
# Voolstra et al PLoS ONE 2011, 6(5): e20392.
library(pheatmap)
# Make sure your input files are in your working directory
setwd("/Users/mikeconnelly/computing/projects/EAPSI_Pocillopora_AxH/R/GO_MWU/")

#Note: must run this script while results of GO_MWU are loaded into R
#must run this to initiate, then can ignore...
gg <- read.table("pdam_emapper_gene2go.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE) #the iso2gene file for P. damicornis

#After running GO_MWU script (download latest version at https://github.com/z0on/GO_MWU)
#rerun the following using the same input/goAnnotations/OBO, etc in the folder where the output files are located (BP_ ; dissim_ etc)

# Edit these as needed, then highlight everything and execute 
input <- "lightcyan1_kme.txt" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="pdam_emapper_gene2go.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision= "BP" # either MF, or BP, or CC
extraOptions='Module=TRUE,Alternative="g"'  # see README_GO_MWU.txt for detail alternative=g for categorical
absValue=-log(0.05,10) # we want to count genes with this or better absolute measure value within each displayed GO category. This does not affect statistics.

level1=0.1 # p-value cutoff for the GO categories to plot as a tree; the worst ones will be printed in small italic font. Specify cutoff=1 to summarize all the genes at or exceeding absValue. 
level2=0.05 # # p-value cutoff for printing in regular font.
level3=0.01 # p-value cutoff for printing in large bold font.
adjusted=T # replace with F to plot un-adjusted p-values.
txtsize=1  # decrease this one to squeeze more GO descriptions on the same panel.
font.family="sans" #"serif"

require(ape)
in.mwu=paste("MWU",goDivision,input,sep="_")
in.dissim=paste("dissim",goDivision,goAnnotations,sep="_")

cutoff=-log(level1,10)
pv=read.table(in.mwu,header=T)
row.names(pv)=pv$term
in.raw=paste(goDivision,input,sep="_")
rsq=read.table(in.raw,sep="\t",header=T)
rsq$term=as.factor(rsq$term)

if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
heat=data.frame(cbind("pval"=pvals)) 
row.names(heat)=pv$term
heat$pval=-log(heat$pval+1e-15,10)
heat$direction=0
heat$direction[pv$delta.rank>0]=1
if (cutoff>0) { 
	goods=subset(heat,pval>=cutoff) 
} else {
	goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
	goods=heat[row.names(heat) %in% goods.names,]
}

colors=c("dodgerblue2","firebrick1","skyblue","lightcoral")
if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) { 
	colors=c("black","black","grey50","grey50")
}
goods.names=row.names(goods)

# reading and subsetting dissimilarity matrix
diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
row.names(diss)=names(diss)
diss.goods = diss[goods.names,goods.names]

# how many genes out of what we started with we account for with our best categories?
good.len=c();good.genes=c()
for (g in goods.names) {
	sel=rsq[rsq$term==g,]	
	pass=abs(sel$value)>=absValue
	sel=sel[pass,]
	good.genes=append(good.genes,as.character(sel$seq))
	good.len=append(good.len,nrow(sel))
}
ngenes=length(unique(good.genes))

#hist(rsq$value)
totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")

# clustering terms better than cutoff
GO.categories=as.dist(diss.goods)
cl.goods=hclust(GO.categories,method="average")
labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
goods=goods[labs,]


goods$pval2<-10^(-(goods$pval))
goods
################################################################################

###############finding what significant GO genes are 

#must change these each time
golist=read.table("BP_lightcyan1_kme.txt",sep="	",header=T) #read in proper GO list from gomwu output - BP/MF/CC and module name
d5=read.csv("host_vsd_MMlightcyan.csv") #read in proper VSD file for genes in module - generated in wgcna script

mat_go  <- assay(vsd)[golist$seq, ]

rownames(mat_go) # check that rownames are gene names
edata=c(2:45) #only columns with your expression data 2:45 host; 2:58 sym

gene=subset(pv, name=="DNA repair") #write GO term of interest here from your sig list#
t=rownames(gene)
t

is=golist$seq[grep(t,golist$term,ignore.case=T)]
length(is)
#should be the same as in the goods table for your term of interest

#####################first loop through genes matching GO term in module or in "good expression" subset
sel=c();gnms=c()
for ( i in is){
	if (i %in% rownames(mat_go)){
		sel=rbind(sel,mat_go[rownames(mat_go)==i,])
		gnms=append(gnms,substr(paste(as.character(gg$V2[gg$V1==i]),collapse="."),1,200))
	}
}
row.names(sel)=paste(gnms,sel$X,sep=".")
nrow(sel)
rownames(sel)


exp=sel
if (length(exp[,1])==1) { exp=rbind(exp,exp) }
nrow(exp)
#should match 'good genes' tabulation in figure
rownames(exp)

head(exp)

write.csv(exp,file="Host_annotDNArepairGenes_lightcyan1.csv",quote=FALSE)

pheatmap(exp)
