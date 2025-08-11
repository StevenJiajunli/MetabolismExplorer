<img width="1536" height="1024" alt="image" src="https://github.com/user-attachments/assets/d8f250fe-d049-4c57-a3ba-13c714ccf195" />


# 🌿 **Welcome to MetabolismExplorer!**

📦 **MetabolismExplorer** is a powerful tool designed for high-throughput metabolic state analysis.

🔍 **We have curated**:  
- **134** 🧬 *metabolism-related gene sets*  
- **3,673** 💥 *metabolite-associated gene sets*

🧪 **What you need:**  
Just provide your **expression matrix** and **sample classification**.

🚀 **What you get:**  
This tool will **automatically identify differential metabolic states**, helping you **uncover key biological insights** and **guide your downstream research** with ease and precision.

<img width="974" height="358" alt="image" src="https://github.com/user-attachments/assets/05b5d488-058d-4e2f-b8d0-f926fcc216e8" />

✅ **Before you start:**  
Before using MetabolismExplorer, **make sure the following R packages are installed**. You can use the following code to automatically install all required dependencies.


```r
# CRAN packages
cran_packages <- c(
  "tibble", "forcats", "ggplot2", "ggrepel", "ggthemes", 
  "gridExtra", "dplyr", "reshape2"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "Pi", "limma", "metaMA", "statmod", 
  "clusterProfiler", "enrichplot", "msigdbr"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}
```


## 📥 Input Format

### 1. Expression Matrix (`input`)

- **Rows**: Genes  
- **Columns**: Sample IDs  
- **Values**: Normalized expression values (e.g., TPM, FPKM, or log2-transformed counts)  
- The row names should be gene symbols, and the column names should match the sample IDs in the `info` table.

Example:

|         | GSM1594219 | GSM1594220 | GSM1594221 |
|---------|------------|------------|------------|
| EEF1A1  | 13.78      | 13.61      | 14.04      |
| GAPDH   | 11.08      | 10.18      | 10.34      |
| SLC35E2 | 3.32       | 3.32       | 3.32       |
| DUSP22  | 3.64       | 3.32       | 3.32       |
| RPS28   | 11.44      | 10.98      | 12.81      |

---

### 2. Sample Group Information (`info`)

A two-column data frame:

- **id**: Sample IDs (should match input column names)  
- **type**: Grouping information (e.g., HC vs. SLE)

Example:

| id         | type |
|------------|------|
| GSM1594219 | HC   |
| GSM1594220 | HC   |
| GSM1594221 | HC   |
| GSM1594222 | SLE  |
| GSM1594223 | SLE  |
| GSM1594224 | SLE  |

---


## 🚀 Usage 1: Exploration of varied metabolic states between two groups

```r
# Load the package
library(MetabolismExplorer)

# Perform differential metabolic analysis
diff <- me_diffprep(
  data = input,
  info = info,
  group1 = "SLE",
  group2 = "HC"
)

# Perform metabolism state analysis
meta <- me_analysis(
  input = diff,
  geneset = "metabolite",  # options: subsystem, metabolite
  set.min = 5,
  set.max = 1000
)

# Volcano plot
me_volcano_plot(input = meta,
                output = "metabolic enrichment.pdf",
                thres.p = 0.05,
                thres.nes = 0.5,
                topn = 4,
                marker = NULL,
                label.size = 5,
                width = 6, height = 6)

# Lollipop plot
lollipop_plot <- me_lollipop_plot(input = meta)

```
---

## 📂 Available Gene Sets

- **subsystem**: 134 metabolism-related gene sets
- **metabolite**: 3,673 metabolite-associated gene sets

---

## 📊 Example Output

**meta**

| id                                     | ES        | NES       | pvalue     | FDR       | setsize | leading_prop | leading_gene                                            |
|---------------------------------------|-----------|-----------|------------|-----------|---------|---------------|---------------------------------------------------------|
| Glycolysis / Gluconeogenesis          | 0.7914154 | 1.4025731 | 0.04935120 | 0.7440891 | 76      | 29%           | GAPDH/ENO1/LDHA/PGK1/ALDOA/ALDOC/...                   |
| Oxidative phosphorylation             | 0.8006821 | 1.3971383 | 0.05361305 | 0.7440891 | 67      | 52%           | COX5A/UQCRFS1/COX5B/COX6C/NDUFB8/...                   |
| Arachidonic acid metabolism           | 0.7583923 | 1.3869413 | 0.05018307 | 0.7440891 | 105     | 13%           | PPT1/GPX1/GSTO1/ALOX5/CYP2J2/...                       |
| Miscellaneous                         | 0.7643188 | 1.3761349 | 0.06001727 | 0.7440891 | 88      | 19%           | CA1/MPO/TXN/CAT/ALPL/TPP1/...                          |
| Pentose phosphate pathway             | 0.8594855 | 1.3708231 | 0.07773196 | 0.7440891 | 29      | 38%           | SAT1/TALDO1/TKT/PGD/DERA/...                          |
| Glutathione metabolism                | 0.8520904 | 1.3637277 | 0.08004952 | 0.7440891 | 30      | 27%           | GPX1/LAP3/PRDX1/PRDX5/GR/...                          |
| Linoleate metabolism                  | 0.7885787 | 1.3594780 | 0.0750371  | 0.7440891 | 61      | 5%            | GPX1/CYP1B1/GPX4                                      |
| Amino sugar and nucleotide sugar metab| 0.8459103 | 1.3593352 | 0.08415435 | 0.7440891 | 31      | 35%           | HEXB/NAGK/HK3/NPL/RENBP/...                           |
| Starch and sucrose metabolism         | 0.8482701 | 1.3576135 | 0.08438209 | 0.7440891 | 30      | 20%           | AMY1C/PYGL/GYG1/GBE1/GAA/...                          |
| Estrogen metabolism                   | 0.7796301 | 1.3574625 | 0.07587076 | 0.7440891 | 66      | 20%           | MPO/GSTO1/CYP1B1/ALDH1A1/GST3/...                     |

<div align="center">

| Volcano plot | Lollipop plot |
|--------------|---------------|
| <img src="https://github.com/user-attachments/assets/4c1e0f1c-17a6-4944-ab0e-433435e830ff" width="400"/> | <img src="https://github.com/user-attachments/assets/569a701c-74dc-4937-8ff3-40f8bd9824de" width="400"/> |

</div>



---


## 🚀 Usage 2: Leading edge analysis

```r
### Note: Please keep Internet connected while running leading edge analysis ~

# Load the package
library(MetabolismExplorer)

# Load gene sets
path <- system.file("extdata", "all_genesets_gem.rds", package = "MetabolismExplorer")
all_genesets <- readRDS(path)

# Perform differential metabolic analysis
diff <- me_diffprep(
  data = input,
  info = info,
  group1 = "SLE",
  group2 = "HC"
)

# Perform leading edge analysis
list <- all_genesets[["subsystem"]]["Pentose and glucuronate interconversions"]
leading <- me_leading_edge(input = diff,
                           geneset = list)

# Visualization
list <- all_genesets[["subsystem"]]["Pentose and glucuronate interconversions"]
select <- leading$id # Or select genes that you interested in
leading_plot <- me_leading_edge_plot(input = diff,
                                     geneset = list,
                                     leading = TRUE,
                                     select.gene = TRUE,
                                     select.name = select)
leading_plot

```

---

## 📊 Example Output

**Leading genes**
| id     | rank  |
|--------|-------|
| TBX6   | 23740 |
| ALG5   | 23741 |
| SORD   | 23918 |
| MECR   | 24394 |
| XYLB   | 25366 |
| UXS1   | 25580 |
| CRYL1  | 27376 |
| DCXR   | 27789 |
| AKR7A2 | 28687 |
| AKR1B1 | 28774 |

**Leading-edge analysis visualization**

<img width="1174" height="812" alt="image" src="https://github.com/user-attachments/assets/bfbbfaf5-4082-4093-9d5c-4e768e96e195" />

---

## 🚀 Usage 3: ssGSEA

To assess pathway-specific activity across individual samples, **single-sample GSEA (ssGSEA)** can be applied using the `GSVA` package.

### 🔧 Required Input
- `input`: A **gene expression matrix** (rows: genes, columns: samples)
- `inst/extdata/all_genesets_gem.rds`: Precompiled gene sets (`subsystem`, `metabolite`)

---

### 🧬 Load Gene Sets and ssGSEA pipeline

```r
# Load gene sets from the built-in RDS file
rds_path <- system.file("extdata", "all_genesets_gem.rds", package = "MetabolismExplorer")
gs <- readRDS(rds_path)

library(GSVA)

# Run ssGSEA
ssgsea_score <- gsva(
  expr = as.matrix(input),
  gset.idx.list = gs[["subsystem"]]["Glutathione metabolism"],
  method = "ssgsea",
  ssgsea.norm = TRUE,
  verbose = TRUE
)
```

