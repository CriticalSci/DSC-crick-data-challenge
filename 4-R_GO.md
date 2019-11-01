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

## GO Analysis with clusterProfiler

```python
library(clusterProfiler)
library("org.Hs.eg.db")
```

Re-load files if necessary:

```python
rna_diffx <- read.table(file="output/rna_diffx.bed", sep = "\t",col.names = c("chr","fstart","fend","name","score","strand"))
rna_cts <- as.matrix(read.csv("data_challenge/rnaseq_gene_counts.txt",sep="\t",row.names="featureid"))
rna_up <- read.table(file="output/rna_up.bed", sep = "\t", col.names = c("chr","fstart","fend","name","score","strand"))
rna_dw <- read.table(file="output/rna_dw.bed", sep = "\t", col.names = c("chr","fstart","fend","name","score","strand"))
head(rna_diffx,2)
head(rna_cts,2)
```

### GO Analysis for Up/Down Regulated Gene Sets

```python
rna_bgd <- row.names(rna_cts[rowSums(rna_cts) > 0,])
rna_up_list <- as.vector(rna_up$name)
rna_dw_list <- as.vector(rna_dw$name)
rna_up_ego <- enrichGO(gene = rna_up_list, universe = rna_bgd, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
rna_dw_ego <- enrichGO(gene = rna_dw_list, universe = rna_bgd, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
```

```python
rna_up_ego_gn <- setReadable(rna_up_ego, 'org.Hs.eg.db', 'ENSEMBL')
rna_dw_ego_gn <- setReadable(rna_dw_ego, 'org.Hs.eg.db', 'ENSEMBL')
```

```python
goplot(rna_up_ego_gn)
goplot(rna_dw_ego_gn)
png("output/plot/rna_up_goplot.png", width = 1600, height = 1600)
goplot(rna_up_ego_gn)
dev.off()
png("output/plot/rna_dw_goplot.png", width = 1600, height = 1600)
goplot(rna_dw_ego_gn)
dev.off()
```

```python
cnetplot(rna_up_ego_gn, circular=TRUE, colorEdge=TRUE, showCategory = 20)
cnetplot(rna_dw_ego_gn, circular=TRUE, colorEdge=TRUE, showCategory = 20)
png("output/plot/rna_up_cnet.png", width = 1600, height = 1600)
cnetplot(rna_up_ego_gn, circular=TRUE, colorEdge=TRUE, showCategory = 20)
dev.off()
png("output/plot/rna_dw_cnet.png", width = 1600, height = 1600)
cnetplot(rna_dw_ego_gn, circular=TRUE, colorEdge=TRUE, showCategory = 20)
dev.off()
```

```python
emapplot(rna_up_ego_gn, showCategory = 100)
emapplot(rna_dw_ego_gn, showCategory = 100)
heatplot(rna_up_ego_gn, showCategory = 50)
heatplot(rna_dw_ego_gn, showCategory = 50)
png("output/plot/rna_up_emap.png", width = 1000, height = 1000)
emapplot(rna_up_ego_gn, showCategory = 100)
dev.off()
png("output/plot/rna_dw_emap.png", width = 1000, height = 1000)
emapplot(rna_dw_ego_gn, showCategory = 100)
dev.off()
png("output/plot/rna_up_heat.png", width = 1600, height = 800)
heatplot(rna_up_ego_gn, showCategory = 50)
dev.off()
png("output/plot/rna_dw_heat.png", width = 1600, height = 800)
heatplot(rna_dw_ego_gn, showCategory = 50)
dev.off()
```

```python
rna_dx_list <- rna_diffx[,5]
names(rna_dx_list) = as.character(rna_diffx[,4])
rna_dx_list = sort(rna_dx_list, decreasing = TRUE)
head(rna_dx_list)
```

```python
rna_gsea <- gseGO(geneList = rna_dx_list, OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL')
```

```python
rna_gsea_short <- setReadable(rna_gsea, 'org.Hs.eg.db', 'ENSEMBL')
head(rna_gsea_short,2)
```

```python
png("output/plot/rna_diffx_gsea.png", width = 9000, height = 800)
heatplot(rna_gsea_short, foldChange = rna_dx_list)
dev.off()
```

### GO Analysis for Up/Down Gene Sets in combination with Hi/Lo Accessibility in
their regulatory regions (-10000 to +1000)

```python
up_hi <- read.table(file="output/upregby_hiprom.bed", sep = "\t",col.names = c("chr","fstart","fend","name","score","strand"))
dw_hi <- read.table(file="output/dwregby_hiprom.bed", sep = "\t",col.names = c("chr","fstart","fend","name","score","strand"))
up_lo <- read.table(file="output/upregby_loprom.bed", sep = "\t",col.names = c("chr","fstart","fend","name","score","strand"))
dw_lo <- read.table(file="output/dwregby_loprom.bed", sep = "\t",col.names = c("chr","fstart","fend","name","score","strand"))
head(dw_hi,2)
```

```python
up_hi_list <- as.vector(up_hi$name)
dw_hi_list <- as.vector(dw_hi$name)
up_lo_list <- as.vector(up_lo$name)
dw_lo_list <- as.vector(dw_lo$name)
up_hi_ego <- enrichGO(gene = up_hi_list, universe = rna_bgd, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
dw_hi_ego <- enrichGO(gene = dw_hi_list, universe = rna_bgd, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
up_lo_ego <- enrichGO(gene = up_lo_list, universe = rna_bgd, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
dw_lo_ego <- enrichGO(gene = dw_lo_list, universe = rna_bgd, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
```

```python
up_hi_ego_gn <- setReadable(up_hi_ego, 'org.Hs.eg.db', 'ENSEMBL')
dw_hi_ego_gn <- setReadable(dw_hi_ego, 'org.Hs.eg.db', 'ENSEMBL')
up_lo_ego_gn <- setReadable(up_lo_ego, 'org.Hs.eg.db', 'ENSEMBL')
dw_lo_ego_gn <- setReadable(dw_lo_ego, 'org.Hs.eg.db', 'ENSEMBL')
```

```python
goplot(up_hi_ego_gn)
goplot(up_lo_ego_gn)
#no data?
#goplot(dw_hi_ego_gn)
#goplot(dw_lo_ego_gn)
png("output/plot/up_hi_goplot.png", width = 900, height = 900)
goplot(up_hi_ego_gn)
dev.off()
png("output/plot/up_lo_goplot.png", width = 900, height = 900)
goplot(up_lo_ego_gn)
dev.off()
```

```python
cnetplot(up_hi_ego_gn, circular=TRUE, colorEdge=TRUE, showCategory = 20)
cnetplot(up_lo_ego_gn, circular=TRUE, colorEdge=TRUE, showCategory = 20)
png("output/plot/up_hi_cnet.png", width = 1600, height = 1600)
cnetplot(up_hi_ego_gn, circular=TRUE, colorEdge=TRUE, showCategory = 20)
dev.off()
png("output/plot/up_lo_cnet.png", width = 1600, height = 1600)
cnetplot(up_lo_ego_gn, circular=TRUE, colorEdge=TRUE, showCategory = 20)
dev.off()
```

```python
heatplot(up_hi_ego_gn, showCategory = 50)
heatplot(up_lo_ego_gn, showCategory = 50)
png("output/plot/up_hi_heat.png", width = 1600, height = 800)
heatplot(up_hi_ego_gn, showCategory = 50)
dev.off()
png("output/plot/up_lo_heat.png", width = 1600, height = 800)
heatplot(up_lo_ego_gn, showCategory = 50)
dev.off()
```

### GO Analysis on Shortlists

```python
stock_anno <- read.csv("data_challenge/rnaseq_annotation.txt",sep="\t",header = TRUE)
colnames(stock_anno)
```

```python
id2name <- stock_anno[c("gene_id","gene_name")]
head(id2name)
```

```python
rna_diffx_anno <- merge(id2name, rna_diffx, by.x = "gene_id", by.y = "name", all.y = TRUE)
```

```python
uphiprox_TFs <- read.csv("output/uphiprox_TFMs.tsv",sep="\t",header = TRUE)
uphiprox_TFs_diffx <- merge(rna_diffx_anno, uphiprox_TFs, by.x = "gene_name", by.y = "motif_alt_ID", all.y = TRUE)
uphiprox_TFs_diffx
```

```python
hi_TFs <- read.csv("output/hi_TFMs.tsv",sep="\t",header = TRUE)
hi_TFs_diffx <- merge(rna_diffx_anno, hi_TFs, by.x = "gene_name", by.y = "motif_alt_ID", all.y = TRUE)
hi_TFs_diffx <- hi_TFs_diffx[!duplicated(hi_TFs_diffx$gene_id),]
hi_TFs_diffx <- na.omit(hi_TFs_diffx)
head(hi_TFs_diffx)
```

```python
lo_TFs <- read.csv("output/lo_TFMs.tsv",sep="\t",header = TRUE)
lo_TFs_diffx <- merge(rna_diffx_anno, lo_TFs, by.x = "gene_name", by.y = "motif_alt_ID", all.y = TRUE)
lo_TFs_diffx <- hi_TFs_diffx[!duplicated(lo_TFs_diffx$gene_id),]
lo_TFs_diffx <- na.omit(lo_TFs_diffx)
lo_TFs_vshi <- merge(lo_TFs_diffx, hi_TFs_diffx, by.x = "gene_id", by.y = "gene_id", all.x = TRUE)
lo_TFs_vshi <- lo_TFs_diffx[is.na(lo_TFs_diffx$gene_name.y),]
lo_TFs_vshi
```

None of the lo_TF sites are unique versus hi_TF sites

```python
hi_TFs_vslo <- merge(hi_TFs_diffx, lo_TFs_diffx, by.x = "gene_id", by.y = "gene_id", all.x = TRUE)
hi_TFs_vslo <- hi_TFs_vslo[is.na(hi_TFs_vslo$gene_name.y),]
```

```python
hi_TFs_vslo_list <- hi_TFs_vslo[,6]
names(hi_TFs_vslo_list) = as.character(hi_TFs_vslo[,1])
hi_TFs_vslo_list = sort(hi_TFs_vslo_list, decreasing = TRUE)
hi_TFs_vslo_vec <- names(hi_TFs_vslo_list)
hi_TFs_vslo_ego <- enrichGO(gene = hi_TFs_vslo_vec, universe = rna_bgd, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
hi_TFs_vslo_ego <- setReadable(hi_TFs_vslo_ego, 'org.Hs.eg.db', 'ENSEMBL')
```

```python
heatplot(hi_TFs_vslo_ego, foldChange=hi_TFs_vslo_list)
cnetplot(hi_TFs_vslo_ego, foldChange=hi_TFs_vslo_list)
png("output/plot/hiTFs_cnet.png", width = 900, height = 900)
cnetplot(hi_TFs_vslo_ego, foldChange=hi_TFs_vslo_list)
dev.off()
png("output/plot/hiTFs_heat.png", width = 900, height = 900)
heatplot(hi_TFs_vslo_ego, foldChange=hi_TFs_vslo_list)
dev.off()
```

## Interpretation
The GO Analysis on the RNA-Seq Data clearly showed many genes
involved in histone methyltransferase activity are downregulated, which may be a
clue to the mechanism of chromatin accessibility changes induced.
The
upregulated genes indicated cytokine/TNF receptors were activated leading to a
cascade of changes triggering changes in phosphorylation activity, chromatin and
transcription factor activation.
These points are especially clear in the emap
plots and heatmaps.

There is also some indication of hormone receptor binding
activity although the cytokine/TNF pathway looks most strongly activated via the
chromatin remodelling as evidenced by the GO terms associated with genes
upregulated with increased accessibility in their regulatory regions
(up_hi_goplot). Perhaps the hormone receptor binding is more upstream of the
chromatin remodelling program and potentially the treatment applied?

The Gene
Set Enrichment Analysis Heatmap (rna_diffx_gsea) also indicated a number of
leukocyte activation genes are regulated. Most notably IL2RA is upregulated
indicating activation of T or B cells. TNF cytokine expression is also highly
upregulated so this may be an immune cell line responding to some pro-
inflammatory treatment? EGR1 is also highly expressed and upstream of DNA
demethylation pathways.

It appears that decrease in chromatin accessibility
liberated some genes involved in cell adhesion from repressors (up_lo_goplot).
The motif search showed some specific TF motifs are enriched in the more
accessible regions. Those TFs found to be differentially expressed by cross
referencing with the RNASeq form a fairly small list: CREB1, FOS, FOSB, HOXA5,
HOXA7, JUNB, JUND, MAFF, MSX1, NFATC2, NFATC3, NFKB2, REL, RELA, TLX2, BACH2,
JDP2, POU6F1, POU3F3.

FOS or c-Fos is the most highly expressed of these TFs,
it's also connected to chromatin binding activity. Together with EGR1 it's a top
candidate stimulated by growth factors (hormones) or cytokine stimuli.
