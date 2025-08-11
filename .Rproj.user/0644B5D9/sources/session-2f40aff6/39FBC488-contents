library(Seurat)
library(plyr)
library(scater)
library(stringr)
library(future)
library(ggplot2)
library(cowplot)
library(reshape2)
library(data.table)
library(tidyverse)
library(future.apply)
library(RColorBrewer)
library(SingleCellExperiment)
library(data.table)
library(tibble)
library(harmony)



mycol <- c(brewer.pal(5,"Set1"), brewer.pal(8,"Set2"),
           brewer.pal(11,"Set3"), brewer.pal(12,"Paired"),
           brewer.pal(8,"Accent"), brewer.pal(11,"Spectral"),
           brewer.pal(11,"BrBG"), brewer.pal(11,"PiYG"),
           brewer.pal(11,"PuOr"),brewer.pal(11,"RdBu"))


source("/home/ug1268/tools/mouse_sc/public_function_mouse.R")
source("/home/ug1268/tools/mouse_sc/citeseq_function.R")
source("/home/ug1268/program/ICB/cohorts/autocluster.R")

## 谷胱甘肽代谢

SLE_cohorts_diffexp <- readRDS("~/Steven Lijiajun/Steven/SLE/SLE_cohorts_diffexp.rds")
GSE65391 <- SLE_cohorts_diffexp[["GSE65391"]]

# ssgsea

library(GSVA)

GSE65391_exp <- GSE65391[,-1]
GSE65391_exp <- t(GSE65391_exp)
data_new <- as.matrix(GSE65391_exp)

gs <- readRDS("~/Steven Lijiajun/交大仁济肿瘤所/细胞系富集/all_genesets_gem.rds")

ssgsea_score = gsva(data_new, gs[["subsystem"]]["Glutathione metabolism"], method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)

# dotplot

ssgsea_score <- as.data.frame(t(ssgsea_score))

input <- data.frame(id = rownames(ssgsea_score),
                    value = ssgsea_score$`Glutathione metabolism`,
                    type = GSE65391$type)

plot <- common_dotbox(input = input, 
                      method = "t.test")
plot

# barcode plot

diff_GSE65391 <- readRDS("~/Steven Lijiajun/Steven/SLE/diff_GSE65391.rds")

## 富集分析输入文件
input <- data.frame(id = diff_GSE65391$id,
                    value = diff_GSE65391$es)

input_vol <- data.frame(id = diff_GSE65391$id,
                        logfc = diff_GSE65391$logfc,
                        pvalue = diff_GSE65391$pvalue)

colnames(input_vol) <- c("id","logfc","pvalue")

select <- c("LAP3","GPX1","OPLAH","PRDX5","ANPEP","GCLC","PRDX6")

plot <- common_volcano_plot(input = input_vol,
                            output = NULL,
                            thres.p = 0.05,
                            thres.fc = 0.25,
                            topn = 5,
                            marker = select,
                            label.size = 4,
                            width = 8, height = 7)
plot

## 读取代谢通路与代谢物相关信息
path = "~/Steven Lijiajun/交大仁济肿瘤所/细胞系富集/all_genesets_gem.rds"

# 代谢物信息
geneset1 = "metabolite"

# 代谢通路信息
geneset2 = "subsystem"

# GSEA可视化1
plot <- gsea_barcode(input = input,
                     source = path,
                     geneset = geneset2,
                     select = "Bile acid recycling")
plot


meta <- Combined_cleaned@meta.data
input <- data.frame(type = meta$celltype_major,
                    value = meta$`Glutathione metabolism`)

type_order <- input %>%
  group_by(type) %>%
  summarise(med = mean(value, na.rm = TRUE)) %>%
  arrange(med) %>%
  pull(type)

# 设置新的因子顺序
input$type <- factor(input$type, levels = type_order)

ggplot(data = input,
       aes(x = value, y = type, fill = type)) +
  geom_density_ridges(alpha = 1, 
                      color = 'white',
                      rel_min_height = 0.02, #尾部修剪，数值越大修剪程度越高
                      scale = 1.8, #山脊重叠程度调整，scale = 1时刚好触及基线，数值越大重叠度越高
                      quantile_lines = TRUE, #显示分位数线
                      quantiles = 2 #仅显示中位数线
  ) +scale_x_continuous(limits = c(0,0.25),
                        breaks = seq(0, 0.25, by = 0.1))+ #x轴限制
  theme_classic() +
  theme(legend.position = 'none')

### 单细胞CD14_Mono 谷胱甘肽代谢重编程

Combined_cleaned <- readRDS("~/Steven Lijiajun/Steven/SLE/GSE135779_combined.rds")

# Glutathione metabolism

plot <- featureplot_new(data = Combined_cleaned,
                        reduction = "umap_harmony",
                        pt.size = 1, 
                        color = "blue2red",
                        features = "Glutathione metabolism",
                        raster = NULL,
                        outlier.rm = TRUE)
plot

## 提取CD14_Mono

table(Combined_cleaned$celltype_major)

CD14Mono <- subset(Combined_cleaned, subset = celltype_major == "c07: CD14_Mono")

# seurat pilepine

nfeatures = 2000
ndim = 15
neigh = 50
dist = 0.5
res = 0.2

CD14Mono <- NormalizeData(CD14Mono, scale.factor = 10000,
                          normalization.method = "LogNormalize")

CD14Mono <- FindVariableFeatures(CD14Mono, nfeatures = nfeatures, 
                                 selection.method = "vst")

CD14Mono <- ScaleData(CD14Mono, features = VariableFeatures(CD14Mono))

CD14Mono <- RunPCA(CD14Mono, assay = 'RNA', slot = 'scale.data')

CD14Mono <- RunHarmony(CD14Mono, group.by.vars = "sample",dims.use = 1:50,
                       assay.use = "RNA")

CD14Mono <- FindNeighbors(CD14Mono, k.param = neigh,
                          dims = 1:ndim, reduction = "harmony")

CD14Mono <- FindClusters(CD14Mono, resolution = 0.3, n.iter = 50)

CD14Mono <- RunUMAP(CD14Mono, dims = 1:ndim,
                    n.neighbors = neigh, min.dist = dist, 
                    reduction = "harmony", reduction.name = "umap_harmony")

plot <- dimplot_new(data = CD14Mono,
                    reduction = "umap_harmony",
                    pt.size = 0.0000001, label = F,
                    group.by = c("group"))
plot

# 计算CD14_Mono的Glutathione metabolism

path = "~/Steven Lijiajun/交大仁济肿瘤所/细胞系富集/all_genesets_gem.rds"
geneset1 = "metabolite"
geneset2 = "subsystem"

subsystem_score <- seurat_score(data = CD14Mono,
                                source = path,
                                geneset = geneset2,
                                min.sz = 5)

CD14Mono <- AddMetaData(CD14Mono, subsystem_score)

plot <- featureplot_new(data = CD14Mono,
                        reduction = "umap_harmony",
                        pt.size = 0.000000000000000001, 
                        color = "blue2red",
                        features = "Glutathione metabolism",
                        raster = NULL,
                        outlier.rm = FALSE)
plot

## 小提琴plot

meta <- CD14Mono@meta.data
input <- data.frame(type = meta$disease,
                    value = meta$`Glutathione metabolism`)

library(dplyr)
library(ggplot2)

# 如果你之前把 input$value 全部变成 NA，这里再转换一下
input$value <- as.numeric(input$value)
input$type <- as.character(input$type)

type_order <- input %>%
  group_by(type) %>%
  summarise(med = mean(value, na.rm = TRUE)) %>%
  arrange(med) %>%
  pull(type)

# 显示排序结果
print(type_order)

# 设置新的因子顺序
input$type <- factor(input$type, levels = type_order)

ggplot(input, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.colour = NA) +
  stat_compare_means(method = "t.test") +  # 添加显著性检验
  theme_classic() +
  theme(
    axis.line = element_line(colour = 'black'),
    axis.text = element_text(colour = 'black'),
    axis.title.x = element_blank(),
    legend.position = 'top',
    text = element_text(family = 'sans')
  )


# 直方图比较
# Calculate the 25% and 75% percentiles

dtt <- data.frame(celltype = CD14Mono$celltype_major, 
                  Pyrimidine = CD14Mono$`Glutathione metabolism`)

low_threshold <- quantile(dtt$Pyrimidine, 0.25)  # 25th percentile
high_threshold <- quantile(dtt$Pyrimidine, 0.75)  # 75th percentile

sum(dtt$Pyrimidine > -0.0482 & dtt$Pyrimidine < 0.0582)
sum(dtt$Pyrimidine < -0.0482)
sum(dtt$Pyrimidine > 0.0582)
summary(dtt$Pyrimidine)

##画图
library(ggplot2)
# 假设数据和分位数已经计算好
# 根据上面数值填充
scoring_counts <- data.frame(
  Group = c("LGlutathione", "DTGlutathione", "HGlutathione"),
  Count = c(14568, 29083, 14563),  #！ 改 填实际的细胞数量
  Boundary = c(-0.15, 0.05, 0.25)  # 数字展示的X轴位置 一般不用改 后期调整即可
)

scoring_counts$Xpos <- c(-0.2, 0, 0.2)

data <- data.frame(Scoing = dtt$Pyrimidine)
quantiles <- c(-0.0482, 0.0582) #！ 改 填25% 和 75% 对应的数值
labels <- c("Score < -0.0482", "Score > 0.0582") #！ 改 填25% 和 75% 对应的数值

#skyblue
# 绘图
Fig2d <- ggplot(data, aes(x = Scoing)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.005, fill = "grey", color = "grey", alpha = 0.4) +
  geom_density(color = "magenta", linetype = "dashed", size = 1) +
  geom_vline(xintercept = quantiles, color = "#BB0021FF" , linetype = "dashed", size = 1) +
  #geom_text(aes(x = Boundary, y = 0.6, label = paste(Group, ':', Count, 'cells')), data = scoring_counts, vjust = -0.5, color = "black", size = 3) +
  geom_text(aes(x = Xpos, y = 3, label = paste(Group, '\n', Count, 'cells')), 
            data = scoring_counts, vjust = -0.5, color = "black", size = 4) +
  geom_text(aes(x = -0.1, y = 1.5, label = labels[1]), hjust = 1, color = "#BB0021FF" , size = 4) +
  geom_text(aes(x = 0.125, y = 1.5, label = labels[2]), hjust = 0, color = "#BB0021FF" , size = 4) +
  labs(title = " ", x = "Scoring", y = "Density") +
  theme_classic()+
  scale_x_continuous(limits = c(-0.3, 0.4))
Fig2d


# 经典红灰绿

summary(CD14Mono$`Glutathione metabolism`)

CD14Mono$Pyri_state <- ifelse(
  CD14Mono$`Glutathione metabolism`< -0.0482, "c01: LGSH",
  ifelse(CD14Mono$`Glutathione metabolism` > 0.0582, "c03: HGSH", "c02: DTGSH")
)

table(CD14Mono$Pyri_state)

plot <- dimplot_new(data = CD14Mono,
                    reduction = "umap_harmony",
                    pt.size = 0.000000000000000001, label = F,
                    group.by = c("Pyri_state"))
plot

# 背对背柱状图 and 棒棒糖图
prop_back2back(datafilt = CD14Mono,
               group = "disease",
               cluster = "Pyri_state",
               order = TRUE)

prop_back2back_lollipop(datafilt = CD14Mono,
                        group = "disease",
                        group1 = "SLE",
                        group2 = "HC",
                        cluster = "Pyri_state")

## Fig2f HGSH vs LGSH

diff_HL <- seurat_diff2(datafilt = CD14Mono,
                        group.by = "Pyri_state",
                        group1 = "c03: HGSH",
                        group2 = "c01: LGSH",
                        assay = "RNA",
                        min.pct = 0.1,
                        thres.fc = 0)

input <- diff_HL[,c(1,2,7)]

colnames(input) <- c("id","logfc","pvalue")

plot <- common_volcano_plot(input = input,
                            output = NULL,
                            thres.p = 0.05,
                            thres.fc = 0.25,
                            topn = 5,
                            marker = "DTYMK",
                            label.size = 4,
                            width = 8, height = 7)
plot

Fig2f <- common_volcano_plot(input = input,
                             output = NULL,
                             thres.p = 0.05,
                             thres.fc = 0.25,
                             topn = 5,
                             marker = NULL,
                             label.size = 4,
                             width = 8, height = 7)
Fig2f

path <- "~/Steven Lijiajun/Steven/to 亲爱的师兄师姐师弟师妹/cqw/Fig.2 嘧啶代谢与DTYMK的引出/all_genesets.rds"

input <- data.frame(id = diff_HL$gene,
                    value =diff_HL$logfc)

result_diffHL_hallmark <- gsea_analysis(input = input,
                                        source = path,
                                        geneset = "msigdb_hallmark",
                                        set.min = 5,
                                        set.max = 1000)

result_diffHL_KEGG <- gsea_analysis(input = input,
                                    source = path,
                                    geneset = "KEGG",
                                    set.min = 5,
                                    set.max = 1000)

result_diffHL_Wiki <- gsea_analysis(input = input,
                                    source = path,
                                    geneset = "GProfiler_Wiki",
                                    set.min = 5,
                                    set.max = 1000)

result_diffHL_KEGG <- result_diffHL_KEGG[order(-result_diffHL_KEGG$NES), ]

rownames(result_diffHL_KEGG) <- 1:nrow(result_diffHL_KEGG)

# lolipop plot

result_diffHL_want <- result_diffHL_KEGG[c(1:4,6,19,286,279,274,283,277),]

plot <- lollipop_plot(result_diffHL_want)

plot

# 选择基因集
list <- all_genesets[["KEGG"]]
list <- list["Systemic lupus erythematosus"]

# leading gene
leading <- leading_edge(input = input, geneset = list)

# 开始可视化
select = c("FCGR1A", "HLA-DRB1", "SSB","SNRPD1")
plot <- leading_edge_plot(input = input,
                          geneset = list,
                          select.gene = TRUE,
                          select.name = select)
plot
