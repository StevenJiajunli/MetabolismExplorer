
#' Title
#'
#' @param input The me_pathwayanalysis result.
#'
#' @return A lollipop plot.
#' @export
#'
#' @examples me_lollipop_plot(input = meta)
me_lollipop_plot <- function(input = meta) {

  library(ggplot2)
  library(forcats)

  # 预处理数据
  input <- data.frame(
    id = input$id,
    NES = input$NES,
    pvalue = input$pvalue
  )

  # 筛选 pvalue < 0.05
  input <- subset(input, pvalue < 0.05)

  # 转换 pvalue 为 -log10
  input$pvalue <- -log10(input$pvalue)

  # 按 NES 排序，并保持因子顺序
  input <- input[order(input$NES),]
  input$id <- fct_inorder(input$id)

  # 开始画图
  ggplot(input, aes(x = id, y = NES, size = pvalue)) +
    geom_segment(aes(x = id, xend = id, y = 0, yend = NES),
                 linewidth = 0.75, color = "grey90") +
    geom_point(aes(color = NES)) +
    geom_hline(yintercept = 0, color = "grey20") +
    coord_flip() +  # 翻转坐标轴
    scale_color_gradient2(low = "#0070b2", high = "#da1735", mid = "white", midpoint = 0) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1)) +
    labs(x = "", y = "Normalized Enrichment Score (NES)", size = "-log10(P-value)")
}
