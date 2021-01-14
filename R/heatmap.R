#' Heatmap of Olink results
#'
#' Given the \code{tibble} object from \code{olink_limma}, plot a heatmap of significant
#' proteins.
#' @param tb A \code{tibble} object from \code{olink_limma} function.
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param p.value \code{numeric}(1): The p.value cutoff for significance level.
#' @param log2FC \code{numeric}(1): The log2FC cutoff for significance level.
#' @param ... graphical parameter passed to \code{pheatmap}.
#' @importFrom pheatmap pheatmap
#' @importFrom scater uniquifyFeatureNames
#' @importFrom magrittr %>%
#' @importFrom grDevices colorRampPalette
#' @export
#' @return A invisible \code{pheatmap} object.
#' @author Ge Tan
#' @examples
#' npxFn <- system.file("extdata", c(
#'   "20200507_Inflammation_NPX_1.xlsx",
#'   "20200625_Inflammation_NPX_2.xlsx"
#' ),
#' package = "OlinkR"
#' )
#' metaFn <- system.file("extdata", "Inflammation_Metadata.xlsx", package = "OlinkR")
#' se <- read_npx(npxFn, metaFn)$SummarizedExperiment
#' tb <- olink_limma(se,
#'   factorCol = "condition_Factor",
#'   contrasts = "Glucose.10mM.Vehicle - Vehicle.Vehicle",
#'   blocking = "Donor_Factor"
#' )
#' olink_heatmap(tb, se)
olink_heatmap <- function(tb, se, p.value = 0.05, log2FC = 0, ...) {
  tb <- tb %>% filter(P.Value < p.value, isPresent, abs(logFC) >= log2FC)
  if (nrow(tb) == 0L) {
    return(NULL)
  }

  toPlot <- tb %>%
    select(intersect(colnames(tb), colnames(se))) %>%
    as.matrix()
  rownames(toPlot) <- uniquifyFeatureNames(tb$OlinkID, tb$Assay)

  annotation_col <- data.frame(colData(se), check.names = FALSE)

  p <- pheatmap(toPlot,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    scale = "row", annotation_col = annotation_col,
    cluster_rows = if_else(nrow(toPlot) >= 2, TRUE, FALSE),
    cluster_cols = FALSE, ...
  )
  invisible(p)
}

#' A overview of NPX value in heatmap
#'
#' Plot the NPX value in a heatmap to spot any potential clustering and outliers
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param scale \code{character}(1): how to scale (z-score) the data. "none": no scaling.
#'                                   "sample": scale across the samples.
#'                                   "feature": scale across the features (proteins).
#' @param cluster_samples \code{logical}(1): do the samples clustering.
#' @param cluster_features \code{logical}(1): do the samples clustering.
#' @param ... graphical parameter passed to \code{pheatmap}.
#' @importFrom pheatmap pheatmap
#' @importMethodsFrom SummarizedExperiment assay
#' @importFrom scater uniquifyFeatureNames
#' @importFrom matrixStats rowSds
#' @importFrom dplyr near
#' @export
#' @return A \code{pheatmap} object.
#' @author Ge Tan
#' @examples
#' npxFn <- system.file("extdata", c(
#'   "20200507_Inflammation_NPX_1.xlsx",
#'   "20200625_Inflammation_NPX_2.xlsx"
#' ),
#' package = "OlinkR"
#' )
#' metaFn <- system.file("extdata", "Inflammation_Metadata.xlsx", package = "OlinkR")
#' se <- read_npx(npxFn, metaFn)$SummarizedExperiment
#' overview_heatmap(se)
overview_heatmap <- function(se, scale = c("none", "row", "column"),
                             cluster_samples = FALSE, cluster_features = FALSE,
                             ...) {
  scale <- match.arg(scale)

  toPlot <- assay(se, "npx")
  rownames(toPlot) <- uniquifyFeatureNames(rownames(se), rowData(se)$Assay)
  annotation_col <- data.frame(colData(se), check.names = FALSE)

  if (scale == "none") {
    p <- pheatmap(toPlot,
      scale = scale, annotation_col = annotation_col,
      cluster_rows = cluster_features, cluster_cols = cluster_samples,
      ...
    )
  } else {
    if (scale == "row") {
      toPlot <- toPlot[!near(rowSds(toPlot), 0), , drop = FALSE]
      p <- pheatmap(toPlot,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        scale = scale, annotation_col = annotation_col,
        cluster_rows = cluster_features, cluster_cols = cluster_samples,
        ...
      )
    } else if (scale == "column") {
      p <- pheatmap(toPlot,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        scale = scale, annotation_col = annotation_col,
        cluster_rows = cluster_features, cluster_cols = cluster_samples,
        ...
      )
    }
  }
  return(p)
}
