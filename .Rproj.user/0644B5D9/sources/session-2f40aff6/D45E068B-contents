

#' Title Leading edge analysis
#'
#' @param input The me_diffprep results.
#' @param geneset The metabolism-related gene set that you interested in.
#'
#' @return The leading edge results.
#' @export
#'
#' @examples me_leading_edge(input = diff, geneset = list)
me_leading_edge <- function(input = diff,
                            geneset = list) {
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
  # 导出结果

  result <- data.frame(rank = eGSEA[["leading"]][[1]])
  result <- rownames_to_column(result, var = "id")
  result
}
