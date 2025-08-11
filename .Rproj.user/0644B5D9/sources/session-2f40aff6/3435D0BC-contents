
#' Title Leading edge analysis visualization.
#'
#' @param input The me_diffprep results.
#' @param geneset The metabolism-related gene set that you interested in.
#' @param leading TRUE
#' @param select.gene TRUE / FALSE
#' @param select.name the metabolism-related genes that you interested in.
#'
#' @return A leading edge plot.
#' @export
#'
#' @examples me_leading_edge_plot(input = diff, geneset = list, leading = TRUE, select.gene = FALSE, select.name = select)
me_leading_edge_plot <- function(input = diff,
                                 geneset = list,
                                 leading = TRUE,
                                 select.gene = FALSE,
                                 select.name = select) {
  library(Pi)
  library(tibble)

  input <- data.frame(id = input$id,
                      value = input$logfc)

  # 数据排序

  input <- input[order(input$value, decreasing = T),]

  # 整理数据

  rownames(input) <- NULL
  input <- column_to_rownames(input, var = "id")
  input$rank <- 1:nrow(input)
  names(input)[1] <- "priority"

  # 开始分析

  eGSEA <- xPierGSEA(input, fast = T,
                     size.range = c(5, 500),
                     nperm = 1000,
                     customised.genesets = geneset)

  # 设定颜色

  col <- "#343391-#343391-#343391-#00b6db-#8dcb8a-#f6bd25-#ea5c2e-#8b2a21-#8b2a21-#8b2a21"
  # col <- "#f3bd2e-#f3bd2e-#f3bd2e-#e2720f-#e02f2c-#861b20-#b51b7f-#7f137f-#622183-#622183-#622183"

  # 开始画图

  if (select.gene == TRUE) {

    plot <- xGSEAdotplot(eGSEA, top = names(geneset),
                         peak.color = "black",
                         leading = leading,
                         leading.query = select.name,
                         leading.query.only = TRUE,
                         leading.size = 3,
                         leading.color = "black",
                         leading.alpha = 1,
                         colormap = col,
                         ncolors = 6)
    plot

  } else {

    plot <- xGSEAdotplot(eGSEA, top = names(geneset),
                         peak.color = "black",
                         leading = leading,
                         leading.query.only = TRUE,
                         leading.size = 3,
                         leading.color = "black",
                         leading.alpha = 1,
                         colormap = col,
                         ncolors = 6)
    plot
  }
}
