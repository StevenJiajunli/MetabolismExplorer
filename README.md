<img width="1536" height="1024" alt="ChatGPT Image 2025年8月7日 17_20_06" src="https://github.com/user-attachments/assets/aa8a6109-2672-4375-8599-67317fa112ef" />







# 🌿 **Welcome to MetabolismExplorer!**

📦 **MetabolismExplorer** is a powerful tool designed for high-throughput metabolic state analysis.

🔍 We have curated:  
- **134** 🧬 *metabolism-related gene sets*  
- **3,673** 💥 *metabolite-associated gene sets*

🧪 **What you need:**  
Just provide your **expression matrix** and **sample classification**.

🚀 **What you get:**  
This tool will **automatically identify differential metabolic states**, helping you **uncover key biological insights** and **guide your downstream research** with ease and precision.

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

## 🚀 Example Usage

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

# Perform pathway enrichment analysis
meta <- me_pathwayanalysis(
  input = diff,
  geneset = "subsystem",  # options: subsystem, metabolism, metabolite
  set.min = 5,
  set.max = 1000
)
```

---

## 📂 Available Gene Sets

- **subsystem**: 134 metabolism-related gene sets
- **metabolite**: 3,673 metabolite-associated gene sets

---

## 📊 Example Output

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

---

## ⚙️ ssGSEA Example

To assess pathway-specific activity across individual samples, **single-sample GSEA (ssGSEA)** can be applied using the `GSVA` package.

### 🔧 Required Input
- `data_new`: A **gene expression matrix** (rows: genes, columns: samples)
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
