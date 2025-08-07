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

| Gene    | GSM1594219 | GSM1594220 | GSM1594221 |
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
  data = diff,
  geneset = "subsystem",  # options: subsystem, metabolism, metabolite
  set.min = 5,
  set.max = 1000
)
```

---

## 📂 Available Gene Sets

- **subsystem**: 134 manually curated metabolism pathways  
- **metabolite**: 3,673 metabolite-associated gene sets  

