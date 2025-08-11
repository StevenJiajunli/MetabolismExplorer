


#' @title Preparing for Metabolic Analysis
#'
#' @param data An expression matrix with genes as row names and samples as column names.
#' @param info A data frame containing sample group information.
#' The first column must be named \code{id}, and it should exactly match the column names (i.e., sample names) of the expression matrix.
#' The second column must be named \code{type}, representing binary classification labels.
#' @param group1 Group1
#' @param group2 Group2
#' @param filter FALSE
#'
#' @return A data frame prepared for downstream metabolic analysis.
#' @export
#'
#' @examples me_diffprep(data = data, info = info, group1 = "R", group2 = "NR", filter = FALSE)
me_diffprep <- function(data = data, 
                       info = info,
                       group1 = "R",
                       group2 = "NR",
                       filter = FALSE) {

  # By default, the comparison is group1 - group2

  library(limma)
  library(metaMA)
  library(statmod)

  # Ensure consistent sample order

  common <- intersect(colnames(data), info$id)
  info <- info[info$id %in% common,]
  data <- data[,info$id]
  data <- as.matrix(data)

  # Remove genes with low variance

  var <- do.call(rbind, lapply(rownames(data), function(i){
    data.frame(id = i, value = var(data[i,]))}))
  data <- data[var$id[var$value > 0],]

  # Build design matrix

  grade <- factor(info[,2], levels = c(group2, group1))
  design <- model.matrix(~0 + grade)

  rownames(design) <- info$id
  colnames(design) <- gsub(pattern = "grade", replacement = "",
                           x = colnames(design))
  # Create contrast matrix

  cont.matrix <- as.matrix(c(-1, 1))
  rownames(cont.matrix) <- c(group2, group1)
  colnames(cont.matrix) <- paste0(group1, " - ", group2)

  # Perform linear fitting

  fit <- lmFit(data, design)
  fit <- contrasts.fit(fit, cont.matrix)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  # Get differentially expressed genes

  diff <- topTable(fit, adjust = 'fdr',
                   coef = 1, n = Inf)

  diff <- data.frame(id = rownames(diff),
                     logfc = diff$logFC,
                     t = diff$t,
                     pvalue = diff$P.Value,
                     FDR = diff$adj.P.Val)

  # Calculate group means

  group_mean <- do.call(rbind, lapply(1:nrow(data), function(i){

    input <- data.frame(id = colnames(data), value = data[i,])
    input <- merge(input, info, by = "id")
    input <- na.omit(input)

    cohen <- (mean(input$value[input$type == group1]) -
                mean(input$value[input$type == group2])) / sd(input$value)

    data.frame(id = rownames(data)[i],
               mean_G1 = mean(input$value[input$type == group1]),
               mean_G2 = mean(input$value[input$type == group2]),
               cohen = cohen)
  }))

  colnames(group_mean)[2:3] <- c(group1, group2)

  # Merge group means with differential expression results

  diff <- merge(group_mean, diff, by = "id")

  # Calculate moderated effect size

  es <- effectsize(fit$t, nrow(info),
                   (fit$df.prior + fit$df.residual))

  es <- data.frame(es = es[,"dprime"],
                   es_var = es[,"vardprime"])
  # Summarize results

  es <- es[diff$id,]
  es <- cbind(diff, es)
  rownames(es) <- NULL
  es

}
