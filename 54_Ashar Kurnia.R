# ==============================================================================
# RNA-Seq differential expression tutorial with DESeq2
# GXA Project: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5244/Results
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(update = TRUE, ask = FALSE)
pacBiocManager::install("DESeq2")
BiocManager::install('EnhancedVolcano')
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)

# ------------------------------------------------------------------------------
# Download data
# ------------------------------------------------------------------------------

# Download raw counts from EBI
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-9479/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(counts)

# Download metadata
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-9479/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(metadata)


# ------------------------------------------------------------------------------
# Wrangle data for DESeq2
# ------------------------------------------------------------------------------


# Counts data wrangling
# DESeq expects the counts to have gene IDs as row names
counts$Gene.ID <- as.character(counts$Gene.ID)
rownames(counts) <- counts$Gene.ID
head(counts)


# Separate gene ID with gene counts
genes = counts[, c("Gene.ID", "Gene.Name"),drop=FALSE]
genes = na.omit(genes) #optional

counts = counts[, -c(1, 2)]
head(counts)


#Metadata Wrangling
# DESeq expects the metadata matrix to have sample IDs in the rownames
metadata$Run <- as.character(metadata$Run)
rownames(metadata) = metadata$Run
head(metadata)

# Only keep columns of interest
metadata = metadata[, c("Sample.Characteristic.stimulus."), drop=FALSE] #Drop are used to make data frame with one column
# Look at metadata to see how the variables change with respect to each other
metadata

# Rename column
colnames(metadata) = "stimulus"


# Turn stimulus/treatment into a factor
metadata$stimulus = factor(metadata$stimulus, levels =c("BMP2", "PBS"))


# ------------------------------------------------------------------------------
# Run DESeq
# ------------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~stimulus)

# Ignore genes with low counts# Ignore genes with low cometadataunts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq
dds <- DESeq(dds)

# Compare expression
res <- results(dds, contrast=c("stimulus", "BMP2", "PBS"), alpha=0.05)
summary(res)

na.omit(res)


# Convert res into dataframe
res_df <- as.data.frame(res)

#Add gene ID column
res_df$Gene.ID <- rownames(res_df)

#Merge with gene name
res_annot <- merge(res_df, genes, by="Gene.ID")
res_annot <- res_annot[
  !is.na(res_annot$Gene.Name) &
    res_annot$Gene.Name != "",
]

#Top 50 gene
res_ordered <- res_annot[order(res_annot$padj), ]
top50 <- head(res_ordered, 50)
top50_genes_df <- as.data.frame(top50)



#
)
# ------------------------------------------------------------------------------
# Visualization
# ------------------------------------------------------------------------------
#Boxplot

# Normalized count from DESeq2 result
norm_counts <- counts(dds, normalized = TRUE)

# Log transform (prevent skewed graph)
log_norm <- log2(norm_counts + 1)

# Make boxplot
# normalized box plot
boxplot(
  log_norm,
  las = 2,
  col = "lightblue",
  main = "Boxplot of Normalized Counts",
  ylab = "Log2 Normalized Counts"
)

# top 10 DEG boxplot
top10_ids <- top50$Gene.ID[1:10]

boxplot(
  log_norm[top10_ids, ],
  las = 2,
  col = "lightgreen",
  main = "Top 10 Differentially Expressed Genes",
  ylab = "Log2 Normalized Counts"
)


# Volcano plot
#Set upregulated and Downregulated gene
up <- res_annot[res_annot$log2FoldChange > 2, ]
down <- res_annot[res_annot$log2FoldChange < 0.5, ]

top_up <- head(up[order(up$padj), ], 25)
top_down <- head(down[order(down$padj), ], 25)

top50 <- rbind(top_up, top_down)

EnhancedVolcano(
  res_annot,
  lab = res_annot$Gene.Name,
  x = 'log2FoldChange',
  y = 'padj',  
  selectLab = top50$Gene.Name,
  pCutoff = 0.05,
  FCcutoff = 1
)


