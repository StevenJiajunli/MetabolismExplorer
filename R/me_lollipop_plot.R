
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

  # Preprocess data
  input <- data.frame(
    id = input$id,
    NES = input$NES,
    pvalue = input$pvalue
  )

  # Filter for pvalue < 0.05
  input <- subset(input, pvalue < 0.05)

  # Transform pvalue to -log10
  input$pvalue <- -log10(input$pvalue)

  # Sort by NES and preserve factor order
  input <- input[order(input$NES),]
  input$id <- fct_inorder(input$id)

  # Plot lollipop chart
  ggplot(input, aes(x = id, y = NES, size = pvalue)) +
    geom_segment(aes(x = id, xend = id, y = 0, yend = NES),
                 linewidth = 0.75, color = "grey90") +
    geom_point(aes(color = NES)) +
    geom_hline(yintercept = 0, color = "grey20") +
    coord_flip() +  # Flip coordinates
    scale_color_gradient2(low = "#0070b2", high = "#da1735", mid = "white", midpoint = 0) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1)) +
    labs(x = "", y = "Normalized Enrichment Score (NES)", size = "-log10(P-value)")
}
