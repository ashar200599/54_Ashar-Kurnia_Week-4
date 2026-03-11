
# Capstone Project of BRSP Omics

![R](https://img.shields.io/badge/R-4.5.2-276DC3?style=flat&logo=r)
![DESeq2](https://img.shields.io/badge/DESeq2-Bioconductor-brightgreen)
![clusterProfiler](https://img.shields.io/badge/clusterProfiler-Bioconductor-yellow)
![DAVID](https://img.shields.io/badge/DAVID-NCBI-pink)
![Dataset](https://img.shields.io/badge/EBI-E--MTAB--9479-pink)


A complete Deseq2 analysis, GO enrichment and KEGG pathway analysis using RNA-seq datasets of OC samples were obtained from EBI (https://www.ebi.ac.uk/gxa/experiments/E-MTAB-9479/Experiment%20Design).  

---

## 🧬 Biological Question

> What 10 top gene expressions changed significantly in BMP signalling in ovarian cancer patients?

| | |
|---|---|
| **Dataset** | E-MTAB-9479 — Homo sapiens,ovarian cancer cell |
| **Comparison** | BMP treatment vs Control |
| **Key finding** | BMP induction caused upregulated ID family gene, which correlated bHLH inhibition leading to ovarian cancer|
| **Paper** |Fukuda et al., Cell Death Discovery 6, no. 1, 1033–1049, 2020 |

---

## Results

### Volcano Plot
![Volcano Plot](result/Volcano Plot.png)

### PCA Plot
![PCA Plot](files/plots/pca_plot.png)

### Heatmap
![Heatmap](files/plots/heatmap_top50_DEGs.png)

---
