


#' @title Metabolic Analysis
#'
#' @param input The me_diffprep results.
#' @param geneset subsystem / metabolite
#' @param set.min 5
#' @param set.max 1000
#'
#' @return A data frame containing the GSEA results, including NES, p-values, FDR, leading-edge genes, etc.
#' @export
#'
#' @examples me_pathwayanalysis(input, geneset = subsystem, set.min = 5, set.max = 1000)
me_analysis <- function(input = diff,
                               geneset = NULL,
                               set.min = 5,
                               set.max = 1000) {

  # Load required packages (recommended to be loaded by the user outside the function)
  library(clusterProfiler)
  library(enrichplot)
  library(gridExtra)
  library(msigdbr)
  library(reshape2)
  options(connectionObserver = NULL)

  # Automatically locate and read the built-in RDS file
  rds_path <- system.file("extdata", "all_genesets_gem.rds", package = "MetabolismExplorer")
  all_genesets <- readRDS(rds_path)

  # Extract the specified gene set
  genesets <- reshape2::melt(all_genesets[[geneset]])
  sig_list <- data.frame(term = genesets[,2], gene = genesets[,1])

  # Prepare expression data
  input <- data.frame(id = input$id,
                      value = input$logfc)

  input <- input[order(input$value, decreasing = TRUE), ]
  data <- input$value
  names(data) <- input$id

  # Perform GSEA enrichment analysis
  kk <- GSEA(data, TERM2GENE = sig_list,
             minGSSize = set.min,
             maxGSSize = set.max,
             pvalueCutoff = 1,
             nPermSimple = 10000,
             eps = 0)

  # Organize results
  result <- data.frame(
    id = kk$Description,
    ES = kk$enrichmentScore,
    NES = kk$NES,
    pvalue = kk$pvalue,
    FDR = kk$p.adjust,
    setsize = kk$setSize,
    leading_prop = kk$leading_edge,
    leading_gene = kk$core_enrichment
  )

  result$leading_prop <- sapply(strsplit(result$leading_prop, ","), function(x) sub(".*=", "", x[1]))

  # Sort by NES and return
  result <- result[order(result$NES, decreasing = TRUE), ]
  return(result)
}
