#' Heatmap of Olink results
#'
#' Given the \code{tibble} object from \code{olink_limma}, plot a heatmap of significant
#' proteins.
#' @param tb A \code{tibble} object from \code{olink_limma} function.
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param p.value \code{numeric}(1): The p.value cutoff for significance level.
#' @param log2FC \code{numeric}(1): The log2FC cutoff for significance level.
#' @param filename \code{character}(1): The filename of heatmap output.
#' @importFrom pheatmap pheatmap
#' @importFrom scater uniquifyFeatureNames
#' @export
#' @return A \code{pheatmap} object.
#' @author Ge Tan
#' @examples
#' npxFn <- system.file("extdata", c("20200507_Inflammation_NPX_1.xlsx",
#'                                   "20200625_Inflammation_NPX_2.xlsx"),
#'                      package="OlinkR")
#' metaFn <- system.file("extdata", "Inflammation_Metadata.xlsx", package="OlinkR")
#' se <- readNPX(npxFn, metaFn)$SummarizedExperiment
#' tb <- olink_limma(se, factorCol="condition_Factor",
#'                   contrasts="Glucose.10mM.Vehicle - Vehicle.Vehicle",
#'                   blocking="Donor_Factor")
#' olink_heatmap(tb, se)

olink_heatmap <- function(tb, se, p.value=0.05, log2FC=0, filename=NA){
  tb <- tb %>% filter(P.Value < p.value, isPresent, abs(logFC) >= log2FC)
  if(nrow(tb) == 0L){
    return(NULL)
  }

  toPlot <- tb %>% select(intersect(colnames(tb), colnames(se))) %>% as.matrix()
  rownames(toPlot) <- uniquifyFeatureNames(tb$OlinkID, tb$Assay)

  annotation_col <- data.frame(colData(se), check.names = FALSE)

  p <- pheatmap(toPlot, color=colorRampPalette(c("blue", "white", "red"))(100),
                scale="row", annotation_col=annotation_col,
                cluster_rows=if_else(nrow(toPlot) >=2, TRUE, FALSE),
                cluster_cols=FALSE, clustering_distance_rows="correlation",
                filename=filename,
                cellwidth=10, cellheight=10)
  return(p)
}
