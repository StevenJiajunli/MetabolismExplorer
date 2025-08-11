
#' Title MetabolismExplorer visualization: Volcano plot
#'
#' @param input The result of me_pathwayanalysis.
#' @param output The name of the output file
#' @param thres.p threshold for p value
#' @param thres.nes thershold for NES
#' @param topn highlight the top n altered metabolic signatures
#' @param marker the metabolic signature you interested in
#' @param label.size label size
#' @param width figure width
#' @param height figure hight
#'
#' @return A volcano plot showing the enrichment results.
#' @export
#'
#' @examples me_volcano_plot(input = meta, output = name, thres.p = 0.05, thres.nes = 1, topn = 3, marker = select, label.size = 5, width = 8, height = 7)

me_volcano_plot <- function(input = meta,
                            output = name,
                            thres.p = 0.05,
                            thres.nes = 1,
                            topn = 3,
                            marker = select,
                            label.size = 5,
                            width = 8, height = 7) {

  library(ggplot2)
  library(ggrepel)
  library(ggthemes)
  library(gridExtra)
  library(dplyr)

  input <- input[,c("id","NES","pvalue")]

  colnames(input) <- c("id","logfc","pvalue")

  # Adjust extremely small P-values
  if (sum(input$pvalue == 0) > 0) {
    input$pvalue[input$pvalue == 0] <- 1e-300
  }

  # Calculate the most differential genes
  input$max <- input$logfc * -log10(input$pvalue)

  # Select significant genes
  sigene_h <- input[input$pvalue < thres.p & input$logfc > thres.nes,]
  sigene_l <- input[input$pvalue < thres.p & input$logfc < -thres.nes,]

  # Select top genes
  toph <- top_n(input, topn, max)
  topl <- top_n(input, -topn, max)
  topgene <- rbind(toph, topl)

  # Select marker genes
  select_gene <- input[input$id %in% marker,]

  # Create plot
  plot <- ggplot(data = input, aes(logfc, -log10(pvalue))) +
    geom_point(alpha = 0.5, size = 3, colour = "grey90") +

    # Add threshold lines
    geom_vline(xintercept = c(thres.nes), color = "black", linetype = "dashed", lwd = 0.75) +
    geom_vline(xintercept = c(-thres.nes), color = "black", linetype = "dashed", lwd = 0.75) +
    geom_hline(yintercept = -log10(thres.p), color = "black", linetype = "dashed", lwd = 0.75) +

    # Significant genes
    geom_point(data = sigene_h, alpha = 0.5, size = 3, color = "#dd897a") +
    geom_point(data = sigene_l, alpha = 0.5, size = 3, color = "#86b6da") +

    # Marker genes
    geom_point(data = select_gene, alpha = 1, size = 3.5, colour = "#fcb040") +

    # Top gene labels
    geom_text_repel(data = topgene, aes(label = id), size = label.size) +

    # Marker gene labels
    geom_text_repel(data = select_gene, aes(label=id), size = label.size) +

    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 15),
          axis.title = element_text(colour = "black", size = 15),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1)) +

    labs(x = 'Log Fold Change', y= 'Log (P-value)', title = '')

  # Output result
  if (is.null(output)) {
    plot
  } else {
    ggsave(output, plot, dpi = 600, width = width, height = height)
  }

}
