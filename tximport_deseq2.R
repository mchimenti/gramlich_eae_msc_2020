## RE-Analysis of Gonsalves/Kuehn mouse MS eye model w/ stem cell treatment 
## Date: 10.22.2018
## Author: Michael Chimenti
## Organism: mm10 / mouse
## Aligners: hisat2 / salmon
## Design: Sham, MS model, MS model + treatment 
## Reps: 6

##########
## Imports
##########

#source("https://bioconductor.org/biocLite.R")
#biocLite("DEGreport")

#negative binomial GLM and related
library('DESeq2')
library('calibrate')
library('tximport')
library('readr')
#annotation
library('biomaRt')
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library('tidyverse')
library('pcaExplorer')
#pathway and gene clusters
library('DEGreport')
#library(pathview)
#library(gage)
#library(gageData)
#library(ggplot2)

#setwd("~/iihg/RNA_seq/kuehn_ankrun_mouse_eye/project_oct2018_RESEQ/") 
setwd("~/collab_proj/kuehn/project_oct2018_RESEQ/")

###########
##Function Defs
###########

get_annotation <- function(dds, biomart_dataset, idtype){
  if(is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
         To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="www.ensembl.org",
                  dataset=biomart_dataset,
                  ensemblRedirect = FALSE)  ## Ensembl biomart redirects to [web mirror]; [[force NO]]
  
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)
  
  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id","gene_name","description")
  
  return(anns2)
}

## Volcano Plot function 
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=ext_gene, cex=textcx, offset=0.3, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

#######################################
## tximport > DESeq2 
#######################################
samples <- read.table("samples.csv", sep=',', header=TRUE)
rownames(samples) <- samples$sample

files <- file.path(getwd(), samples$sname, 'salmon', 'quant.sf')
names(files) <- samples$sample

tx2gene <- read_csv(file.path(getwd(), "tx2gene.csv"), col_names = FALSE)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
ddsTxi <- DESeq(ddsTxi)

anno <- get_annotation(ddsTxi, 'mmusculus_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)

rldTxi <- rlog(ddsTxi, blind=FALSE)
pcaExplorer(dds=ddsTxi,annotation=anno,rlt=rldTxi)

## look at dispersion estimates 
plotDispEsts(ddsTxi)
plotMA(object = ddsTxi, alpha = 0.05)
plotPCA(object = rldTxi, intgroup = 'condition')

#looking at PCs 3 and 4
rld_mat <- assay(rldTxi)
pca <- prcomp(t(rld_mat))
df <- cbind(samples, pca$x)
ggplot(df) + geom_point(aes(x=PC3,y=PC4, color = condition))

# drop samples?  
## sample1 looks like an extreme outlier 
ddsTxi <- ddsTxi[ , ddsTxi$sample != 's1']
ddsTxi$sample <- droplevels(ddsTxi$sample)
ddsTxi <- DESeq(ddsTxi)

##DE testing 


###
res_eae <- lfcShrink(ddsTxi, contrast = c("condition","eae","sham"))# type = 'ashr')
res_eae <- na.omit(res_eae)  #drop NA rows
res_eae_sig <- res_eae[res_eae$padj < 0.1 & res_eae$baseMean > 5.0,]
res_eae_ord <- res_eae_sig[order(res_eae_sig$padj),]
res_eae_ord$ext_gene <- anno[row.names(res_eae_ord), "gene_name"]

png("volcano_eae_DEgenes_fdr_10percent.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_eae_ord, main = "Volcano: DE genes in EAE cells vs. sham", lfcthresh=0.3, 
            sigthresh=0.1, textcx=.6, xlim=c(-1.5, 1.5), ylim = c(4,11))
dev.off()

degPlot(dds = ddsTxi, res = res_eae_ord, n = 9, xs = 'condition')
mycols <- c("baseMean", "log2FoldChange", "padj", "ext_gene")
write.csv(x = res_eae_ord[,mycols], file = "eae_vs_sham_DEgenes_FDR_10percent.csv")

###
res_eae_msc_sham <- lfcShrink(ddsTxi, contrast = c("condition","eae_msc","sham"))#type = 'ashr')
res_eae_msc_sham <- na.omit(res_eae_msc_sham)  #drop NA rows
res_eae_msc_sham_sig <- res_eae_msc_sham[res_eae_msc_sham$padj < 0.05 & res_eae_msc_sham$baseMean > 5.0,]
res_eae_msc_sham_ord <- res_eae_msc_sham_sig[order(res_eae_msc_sham_sig$padj),]
res_eae_msc_sham_ord$ext_gene <- anno[row.names(res_eae_msc_sham_ord), "gene_name"]

png("volcano_eae_msc_vs_sham_DEgenes_FDR_5percent.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_eae_msc_sham_ord, main = "Volcano: DE genes in EAE w/ MSC treatment vs. SHAM", lfcthresh=0.5, sigthresh=0.05, textcx=.4, xlim=c(-2, 2), ylim = c(4,10))
dev.off()

degPlot(dds = ddsTxi, res = res_eae_msc_sham_ord, n = 9, xs = "condition")
write.csv(x = res_eae_msc_sham_ord[,mycols], file = "eae_msc_vs_sham_DE_genes_FDR_5percent.csv")

###
res_eae_msc <- lfcShrink(ddsTxi, contrast = c("condition","eae_msc","eae"), type = 'normal')
res_eae_msc <- na.omit(res_eae_msc)  #drop NA rows
res_eae_msc_sig <- res_eae_msc[res_eae_msc$padj < 0.05 & res_eae_msc$baseMean > 5.0,]
res_eae_msc_ord <- res_eae_msc_sig[order(res_eae_msc_sig$padj),]
res_eae_msc_ord$ext_gene <- anno[row.names(res_eae_msc_ord), "gene_name"]
hist(res_eae_msc_ord$log2FoldChange)

png("volcano_eae_msc_vs_eae_FDR_5percent.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_eae_msc_ord, main = "Volcano: DE Genes btw EAE_MSC and EAE, FDR < 0.05", lfcthresh=0.5, sigthresh=0.05, textcx=.4, xlim=c(-2, 2), ylim = c(3,10))
dev.off()

degPlot(dds = ddsTxi, res = res_eae_msc_ord, n = 9, xs = "condition")
write.csv(x = res_eae_msc_ord[,mycols], file = "eae_msc_vs_eae_DE_genes_FDR5percent.csv")


##### Write files for iPathwayGuide

res_eae_msc_ipath <- results(ddsTxi, contrast = c("condition", "eae_msc", "eae"))
res_eae_msc_ipath <- na.omit(res_eae_msc_ipath)
res_eae_msc_ipath$gene <- anno[row.names(res_eae_msc_ipath), "gene_name"]
write.table(as.data.frame(res_eae_msc_ipath), file = "res_eae_msc_for_iPath.txt", sep = '\t', row.names=FALSE)


res_eae_msc_sham_ipath <- results(ddsTxi, contrast = c("condition", "eae_msc", "sham"))
res_eae_msc_sham_ipath <- na.omit(res_eae_msc_sham_ipath)
res_eae_msc_sham_ipath$gene <- anno[row.names(res_eae_msc_sham_ipath), "gene_name"]
write.table(as.data.frame(res_eae_msc_sham_ipath), file = "res_eae_msc_sham_for_iPath.txt", sep = '\t', row.names=FALSE)
