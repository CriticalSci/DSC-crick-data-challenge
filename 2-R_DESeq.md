# Recruitment Data Challenge

The Bioinformatics & Biostatistics Group @ The
Francis Crick Institute

## Introduction

Here you will find the data from an
RNA-Seq and ATAC-Seq experiment. Both experiments have the same design. There is
a treatment and control group each containing three replicates making a total of
six samples per experiment. The data files are defined as follows (all files are
tab delimited text files):

#### RNA-Seq Data

- **rnaseq_design.txt**: Sample
ids and corresponding condition labels.
- **rnaseq_gene_counts.txt**: Raw (not
normalised) gene-level read counts for each sample.
- **rnaseq_annotation.txt**:
Gene level annotation.

#### ATAC-Seq Data

- **atacseq_design.txt**: Sample ids
and corresponding condition labels.
- **atacseq_peak_counts.txt**: Raw (not
normalised) ATAC-Seq peak level counts for each sample.
- **atacseq_peaks.bed**:
A bed file defining the peak loci

All sequence data were aligned to the human
genome reference hg38.

## The Challenge

The treatment here is thought to
activate a transcriptional program via remodelling of the chromatin
architecture. The aim here is to:  
1. Identify genes that may be regulated in
this fashion.  
2. Identify the possible transcriptional programs involved.  
3.
Present candidate transcription factors that may be responsible for the
underlying regulation.  

Please produce a 20 minute presentation detailing your
exploration of the data, your analysis approach and findings?

# Analysis
## Strategy

1. Identify genes with significant changes in
expression.
2. Identify zones with significant changes in accessibility.
3.
Detect hotspots in accessibility changes over gene regulatory areas of
differentially expressed genes.
4. Detect enriched TF motifs in zones presenting
accessibility changes.
5. Detect enriched TF motifs in hotspots.
6. Perform GO
Analysis to put genes in context.

## Differential Expression with DESeq2
### Setup

```python
library("DESeq2")
library("IHW")
```

import data

```python
#atac_cts <- as.matrix(read.csv("data_challenge/atacseq_peak_counts.txt",sep="\t",row.names="peakid"))
atac_fcts <- as.matrix(read.csv("output/atacseq_fpeak_counts.txt",
                                sep="\t",
                                row.names="peakid"))
atac_col <- read.csv("data_challenge/atacseq_design.txt",
                     sep="\t",
                     row.names="sample.id")
dim(atac_fcts)
head(atac_fcts,2)
atac_col
```

```python
rna_cts <- as.matrix(read.csv("data_challenge/rnaseq_gene_counts.txt",
                              sep="\t",
                              row.names="featureid"))
rna_col <- read.csv("data_challenge/rnaseq_design.txt",
                    sep="\t",
                    row.names="sample.id")
dim(rna_cts)
head(rna_cts,2)
rna_col
```

Double check for entry into DESeq2

```python
all(rownames(atac_col) == colnames(atac_fcts))
all(rownames(rna_col) == colnames(rna_cts))
```

```python
atac_DDS <- DESeqDataSetFromMatrix(countData = atac_fcts, 
                                   colData = atac_col, 
                                   design = ~ condition)
atac_DDS
```

```python
rna_DDS <- DESeqDataSetFromMatrix(countData = rna_cts, 
                                  colData = rna_col, 
                                  design = ~ condition)
rna_DDS
```

Set reference

```python
atac_DDS$condition <- relevel(atac_DDS$condition, 
                              ref = "control")
rna_DDS$condition <- relevel(rna_DDS$condition, 
                             ref = "control")
```

### rLog normalisation and PCA

```python
atac_DDS_Rlog <- rlog(atac_DDS)
rna_DDS_Rlog <- rlog(rna_DDS)
par(mfrow=c(1,2))
plotPCA(atac_DDS_Rlog, intgroup = "condition")
plotPCA(rna_DDS_Rlog, intgroup = "condition")
png("output/plot/atac_pca.png", width = 800, height = 800)
plotPCA(atac_DDS_Rlog, intgroup = "condition")
dev.off()
png("output/plot/rna_pca.png", width = 800, height = 800)
plotPCA(rna_DDS_Rlog, intgroup = "condition")
dev.off()
```

### DESeq2 Analysis

```python
atac_DDS_dea <- DESeq(atac_DDS)
rna_DDS_dea <- DESeq(rna_DDS)
```

Quick look at results with Benjamini-Hochberg (BH) FDR < 0.05

```python
atac_DDS_dea_res <- results(atac_DDS_dea)
summary(atac_DDS_dea_res)
sum(atac_DDS_dea_res$padj < 0.05, na.rm = TRUE)
```

```python
rna_DDS_dea_res <- results(rna_DDS_dea)
summary(rna_DDS_dea_res)
sum(rna_DDS_dea_res$padj < 0.05, na.rm = TRUE)
```

Quick look at results with Independent Hypothesis Weighting (IHW) FDR < 0.05

```python
atac_DDS_dea_resIHW <- results(atac_DDS_dea, filterFun=ihw)
summary(atac_DDS_dea_resIHW)
sum(atac_DDS_dea_resIHW$padj < 0.05, na.rm = TRUE)
```

```python
rna_DDS_dea_resIHW <- results(rna_DDS_dea, filterFun=ihw)
summary(rna_DDS_dea_resIHW)
sum(rna_DDS_dea_resIHW$padj < 0.05, na.rm = TRUE)
```

### Plotting

Shrinkage

```python
#check the coef term
resultsNames(atac_DDS_dea)
resultsNames(rna_DDS_dea)
```

```python
atac_DDS_dea_resLFC <- lfcShrink(atac_DDS_dea, coef = "condition_treated_vs_control", type = "apeglm")
rna_DDS_dea_resLFC <- lfcShrink(rna_DDS_dea, coef = "condition_treated_vs_control", type = "apeglm")
```

Scatter plots

```python
par(mfrow=c(2,2))
plotDispEsts(atac_DDS_dea)
plotDispEsts(rna_DDS_dea)
plotMA(atac_DDS_dea_resLFC, main = "ATAC-Seq", alpha = 0.05)
plotMA(rna_DDS_dea_resLFC, main = "RNA-Seq", alpha = 0.05)
#export as png
png("output/plot/diffx.png", width = 1600, height = 1600)
par(mfrow=c(2,2))
plotDispEsts(atac_DDS_dea)
plotDispEsts(rna_DDS_dea)
plotMA(atac_DDS_dea_resLFC, main = "ATAC-Seq", alpha = 0.05)
plotMA(rna_DDS_dea_resLFC, main = "RNA-Seq", alpha = 0.05)
dev.off()
```

### Merge and export results

Intersect genes significantly differentially
expressed in both tests (BH & IHW)

```python
#atacseq
atac_DDF_BH <- as.data.frame(subset(atac_DDS_dea_res, padj < 0.05))
atac_DDF_IHW <- as.data.frame(subset(atac_DDS_dea_resIHW, padj < 0.05))
atac_DDF_FDR <- merge(atac_DDF_BH, atac_DDF_IHW, by=0)
#rnaseq
rna_DDF_BH <- as.data.frame(subset(rna_DDS_dea_res, padj < 0.05))
rna_DDF_IHW <- as.data.frame(subset(rna_DDS_dea_resIHW, padj < 0.05))
rna_DDF_FDR <- merge(rna_DDF_BH, rna_DDF_IHW, by=0)
head(atac_DDF_FDR,2)
dim(atac_DDF_FDR)
head(rna_DDF_FDR,2)
dim(rna_DDF_FDR)
```

```python
#removing redundant columns
atac_DDF_FDR <- atac_DDF_FDR[-c(4:5,8:12)]
rna_DDF_FDR <- rna_DDF_FDR[-c(4:5,8:12)]
head(atac_DDF_FDR,2)
head(rna_DDF_FDR,2)
```

Merging RNA DiffX Results

```python
#merge with coordinates keeping only significant diffx
rna_bed <- read.csv("output/rnaseq_genes.bed6",
                    sep="\t",
                    header = FALSE, 
                    col.names = c("chr","fstart","fend","name","score","strand"))
rna_diffx <- merge(rna_bed, 
                   rna_DDF_FDR, 
                   by.x = "name", 
                   by.y = "Row.names", 
                   all = FALSE)
rna_diffx <- rna_diffx[c(2:4,1,8,6)] #back to bed6 format
head(rna_bed,2)
head(rna_diffx,2)
write.table(rna_diffx, file="output/rna_diffx.bed", 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
```

Merging ATAC DiffA Results

```python
#merge with coordinates keeping only significant diffx
atac_bed <- read.csv("output/atacseq_fpeaks.bed6",sep="\t",header = FALSE, col.names = c("chr","fstart","fend","name","score","strand"))
atac_diffx <- merge(atac_bed, 
                    atac_DDF_FDR, 
                    by.x = "name", 
                    by.y = "Row.names", 
                    all.x = FALSE)
atac_diffx <- atac_diffx[c(2:4,1,8,6)] #back to bed6 format
dim(atac_bed)
head(atac_bed,2)
dim(atac_diffx)
head(atac_diffx,2)
write.table(atac_diffx, file="output/atac_diffx.bed", 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
```

```python
#split and export
atac_hi <- atac_diffx[atac_diffx$log2FoldChange.x > 0,]
write.table(atac_hi, 
            file="output/atac_hi.bed", 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
atac_lo <- atac_diffx[atac_diffx$log2FoldChange.x < 0,]
write.table(atac_lo, 
            file="output/atac_lo.bed", 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
rna_up <- rna_diffx[rna_diffx$log2FoldChange.x > 0,]
write.table(rna_up, 
            file="output/rna_up.bed", 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
rna_dw <- rna_diffx[rna_diffx$log2FoldChange.x < 0,]
write.table(rna_dw, 
            file="output/rna_dw.bed", 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
```

Cleaning Up

```python
#don't need these anymore
detach("package:IHW", unload=TRUE)
detach("package:DESeq2", unload=TRUE)
```

```python
sort(sapply(ls(),function(x){object.size(get(x))}))
```

## ChIPseeker Analysis / Hotspot Identification

```python
#load libraries
library(biomaRt)
library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
```

Grabbing transcript level info via Biomart

```python
listMarts(host="www.ensembl.org")
```

```python
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
searchDatasets(mart = ensembl, pattern = "GRCh38")
```

```python
ensemblHs38 = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene2transcript <- getBM(mart = ensemblHs38, attributes = c("ensembl_gene_id","ensembl_transcript_id"))
head(gene2transcript)
dim(gene2transcript)
```

We can generate custom TxDb objects for up/down regulated genes

```python
rna_up_tids <- merge(rna_up, gene2transcript, by.x = "name", by.y = "ensembl_gene_id", all.x = TRUE)
rna_up_tids <- na.omit(rna_up_tids[,"ensembl_transcript_id"])
head(rna_up_tids)
```

```python
rna_dw_tids <- merge(rna_dw, gene2transcript, by.x = "name", by.y = "ensembl_gene_id", all.x = TRUE)
rna_dw_tids <- na.omit(rna_dw_tids[,"ensembl_transcript_id"])
head(rna_dw_tids)
```

```python
#this may take a while...
txdb_up <- makeTxDbFromBiomart(biomart = "ensembl", dataset="hsapiens_gene_ensembl", transcript_ids=rna_up_tids)
txdb_dw <- makeTxDbFromBiomart(biomart = "ensembl", dataset="hsapiens_gene_ensembl", transcript_ids=rna_dw_tids)
```

Import the data

```python
atac_hi_GR <- readPeakFile("output/atac_hi.bed", as = "GRanges")
atac_lo_GR <- readPeakFile("output/atac_lo.bed", as = "GRanges")
```

```python
covplot(atac_hi_GR, weightCol = "V5") #nice but not so informative
```

Annotating peaks with information regarding their distance/overlaps with up or
down regulated gene features

```python
atac_hi_upreg <- annotatePeak(atac_hi_GR, 
                              tssRegion = c(-10000,1000), 
                              TxDb = txdb_up, 
                              annoDb = "org.Hs.eg.db", #adds extra info like ENTREZ ID and symbol
                              genomicAnnotationPriority = c("Promoter", 
                                                            "5UTR", 
                                                            "3UTR", 
                                                            "Exon", 
                                                            "Intron", 
                                                            "Downstream"))
```

```python
atac_hi_dwreg <- annotatePeak(atac_hi_GR, 
                              tssRegion = c(-10000,1000), 
                              TxDb = txdb_dw, 
                              annoDb = "org.Hs.eg.db", #adds extra info like ENTREZ ID and symbol
                              genomicAnnotationPriority = c("Promoter", 
                                                            "5UTR", 
                                                            "3UTR", 
                                                            "Exon",
                                                            "Intron", 
                                                            "Downstream"))
```

```python
atac_lo_upreg <- annotatePeak(atac_lo_GR, 
                              tssRegion = c(-10000,1000), 
                              TxDb = txdb_up, 
                              annoDb = "org.Hs.eg.db", #adds extra info like ENTREZ ID and symbol
                              genomicAnnotationPriority = c("Promoter", 
                                                            "5UTR", 
                                                            "3UTR", 
                                                            "Exon",
                                                            "Intron", 
                                                            "Downstream"))
```

```python
atac_lo_dwreg <- annotatePeak(atac_lo_GR, 
                              tssRegion = c(-10000,1000), 
                              TxDb = txdb_dw, 
                              annoDb = "org.Hs.eg.db", #adds extra info like ENTREZ ID and symbol
                              genomicAnnotationPriority = c("Promoter", 
                                                            "5UTR", 
                                                            "3UTR", 
                                                            "Exon",
                                                            "Intron", 
                                                            "Downstream"))
```

We can plot an overview of where the accessibility peaks land (distal/intergenic
peaks have been purposefully left out but these simply appear as NAs)

```python
PeakAnnoList <- list("High Upreg"=atac_hi_upreg, 
                     "High Dwreg"=atac_hi_dwreg, 
                     "Low Upreg"=atac_lo_upreg, 
                     "Low Dwreg"=atac_lo_dwreg)
plotAnnoBar(PeakAnnoList)
#export as png
png("output/plot/atac_prox_bar.png", width = 1600, height = 1600)
plotAnnoBar(PeakAnnoList)
dev.off()
```

These steps allow us to build the coverage profiles in the "promoter" region
which we initially defined as -10000 to +1000 from TSS/GeneStart

```python
upreg_prom_locs <- getPromoters(TxDb = txdb_up, 
                                upstream = 10000, 
                                downstream = 1000)
dwreg_prom_locs <- getPromoters(TxDb = txdb_dw, 
                                upstream = 10000, 
                                downstream = 1000)
atac_hi_upreg_prom <- getTagMatrix(atac_hi_GR, 
                                   windows = upreg_prom_locs)
atac_hi_dwreg_prom <- getTagMatrix(atac_hi_GR, 
                                   windows = dwreg_prom_locs)
atac_lo_upreg_prom <- getTagMatrix(atac_lo_GR, 
                                   windows = upreg_prom_locs)
atac_lo_dwreg_prom <- getTagMatrix(atac_lo_GR, 
                                   windows = dwreg_prom_locs)
PromTagMList <- list("High Upreg"=atac_hi_upreg_prom, 
                     "High Dwreg"=atac_hi_dwreg_prom, 
                     "Low Upreg"=atac_lo_upreg_prom, 
                     "Low Dwreg"=atac_lo_dwreg_prom)
```

Heatmap views

```python
tagHeatmap(PromTagMList, 
           xlim=c(-10000, 1000), 
           color=NULL)
#export as png
png("output/plot/atac_prom_heatmap.png", width = 1600, height = 1600)
tagHeatmap(PromTagMList, xlim=c(-10000, 1000), color=NULL)
dev.off()
```

Next is the coverage plot, to identify hotspots visually

```python
plotAvgProf(PromTagMList, 
            xlim = c(-10000,1000), 
            facet="row", 
            ylab = "Accessibility") + scale_x_continuous(breaks = round(seq(-10000, 1000, by = 1000),1))
#export as png
png("output/plot/atac_prom_graph.png", width = 1600, height = 1600)
plotAvgProf(PromTagMList, xlim = c(-10000,1000), facet="row", ylab = "Accessibility") + scale_x_continuous(breaks = round(seq(-10000, 1000, by = 1000),1))
dev.off()
```

based on this output we can select Regions of Interest:  

**High Access /
Upregulated**
- Window 0: -2000 to 0  
- Window 1: -4250 to -3250  
- Window 2:
-8000 to -4500
- Window 3: -9750 to -8500  

**High Access / Downregulated**
-
Window 0: -2000 to 0  
- Window 1: -9500 to -4500  

**Low Access /
Upregulated**
- Window 0: -3000 to 0  
- Window 1: -6000 to -4500  
- Window 2:
-9750 to -8500

**Low Access / Downregulated**
- Window 0: -3000 to 0   
-
Window 1: -9750 to -8250

```python
atac_hi_upregDF <- as.data.frame(atac_hi_upreg)
atac_hi_dwregDF <- as.data.frame(atac_hi_dwreg)
atac_lo_upregDF <- as.data.frame(atac_lo_upreg)
atac_lo_dwregDF <- as.data.frame(atac_lo_dwreg)
write.table(atac_hi_upregDF, file="output/atac_hi_upreg.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(atac_hi_dwregDF, file="output/atac_hi_dwreg.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(atac_lo_upregDF, file="output/atac_lo_upreg.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(atac_lo_dwregDF, file="output/atac_lo_dwreg.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
```

Extracting the gene and coordinates that had hotspots in the "promoter" region

```python
atac_hi_upreg_promDF <- atac_hi_upregDF[grep("Promoter", atac_hi_upregDF$annotation),]
atac_hi_dwreg_promDF <- atac_hi_dwregDF[grep("Promoter", atac_hi_dwregDF$annotation),]
atac_lo_upreg_promDF <- atac_lo_upregDF[grep("Promoter", atac_lo_upregDF$annotation),]
atac_lo_dwreg_promDF <- atac_lo_dwregDF[grep("Promoter", atac_lo_dwregDF$annotation),]
upregby_hiprom <- merge(rna_diffx, atac_hi_upreg_promDF, by.x = "name", by.y = "geneId", all = FALSE)
dwregby_hiprom <- merge(rna_diffx, atac_hi_dwreg_promDF, by.x = "name", by.y = "geneId", all = FALSE)
upregby_loprom <- merge(rna_diffx, atac_lo_upreg_promDF, by.x = "name", by.y = "geneId", all = FALSE)
dwregby_loprom <- merge(rna_diffx, atac_lo_dwreg_promDF, by.x = "name", by.y = "geneId", all = FALSE)
upregby_hiprom <- upregby_hiprom[!duplicated(upregby_hiprom$name),]
dwregby_hiprom <- dwregby_hiprom[!duplicated(dwregby_hiprom$name),]
upregby_loprom <- upregby_loprom[!duplicated(upregby_loprom$name),]
dwregby_loprom <- dwregby_loprom[!duplicated(dwregby_loprom$name),]
#back to bed6 format
upregby_hiprom <- upregby_hiprom[c(2:4,1,5,6)]
dwregby_hiprom <- dwregby_hiprom[c(2:4,1,5,6)]
upregby_loprom <- upregby_loprom[c(2:4,1,5,6)]
dwregby_loprom <- dwregby_loprom[c(2:4,1,5,6)]
dim(upregby_hiprom)
dim(dwregby_hiprom)
dim(upregby_loprom)
dim(dwregby_loprom)
#head(atac_hi_upreg_promDF,2)
#head(upregby_hiprom,2)
write.table(upregby_hiprom, file="output/upregby_hiprom.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(dwregby_hiprom, file="output/dwregby_hiprom.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(upregby_loprom, file="output/upregby_loprom.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(dwregby_loprom, file="output/dwregby_loprom.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
```
