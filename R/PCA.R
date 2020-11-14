#' PCA analysis of olink data
#'
#' Given a \code{SummarizedExperiment} object, perform the PCA analysis implemented
#' in Bioconductor package \code{PCAtools}
#'
#' @param se An \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param colby \code{character}(1): The metadata column in \code{se} to use for colour.
#' @param shape \code{character}(1): The metadata column in \code{se} to use for shape.
#' @param metavars \code{character}(n): The metedata columns with numeric values in \code{se} to correlate with principal components.
#' @param removeVar \code{numeric}(1): Remove this fraction of variables based on low variance. DEFAULT = 0.1. OPTIONAL.
#' @importFrom PCAtools pca parallelPCA findElbowPoint screeplot pairsplot biplot plotloadings eigencorplot getComponents
#' @importFrom scater uniquifyFeatureNames
#' @importFrom ggplot2 geom_text aes
#' @importMethodsFrom SummarizedExperiment assay colData
#' @importFrom utils head
#' @export
#' @return A list of \code{ggplot} and \code{tibble} objects.
#' @author Ge Tan
#' @examples
#' npxFn <- system.file("extdata", c("20200507_Inflammation_NPX_1.xlsx",
#'                                   "20200625_Inflammation_NPX_2.xlsx"),
#'                      package="OlinkR")
#' metaFn <- system.file("extdata", "Inflammation_Metadata.xlsx", package="OlinkR")
#' se <- read_npx(npxFn, metaFn)$SummarizedExperiment
#' ans <- olink_pca(se, colby="condition_Factor", shape="Donor_Factor", metavars="Weight_Numeric")

olink_pca <- function(se, colby=NULL, shape=NULL, metavars=NULL,
                      removeVar=0.1){
  x <- assay(se, "npx")
  rownames(x) <- uniquifyFeatureNames(rownames(se), rowData(se)$Assay)
  while(class(p <- try(pca(x, metadata = colData(se), scale=TRUE,
                           removeVar = removeVar), silent=TRUE)) == "try-error"){
    removeVar <- removeVar + 0.1
  }
  ans <- list()
  ## scree plot
  horn <- parallelPCA(x)
  elbow <- findElbowPoint(p$variance)
  ans[["scree plot"]] <- screeplot(p, components = getComponents(p, 1:20),
                                   vline = c(horn$n, elbow)) +
    geom_text(aes(horn$n + 1, 50, label = "Horn's", vjust = -1)) +
    geom_text(aes(elbow + 1, 50, label = "Elbow", vjust = -1))

  ## pairs plot
  ans[["pairs plot"]] <- pairsplot(p, components = getComponents(p, c(1:3)),
                                   triangle = TRUE, trianglelabSize = 12,
                                   hline = 0, vline = 0,
                                   pointSize = 1, gridlines.major = FALSE,
                                   gridlines.minor = FALSE,
                                   colby = colby, shape = shape,
                                   title="Pairs plot", titleLabSize = 16,
                                   plotaxes = FALSE, returnPlot = TRUE)

  ## bi-plots
  ans[["bi-plot"]] <- biplot(p, colby = colby, shape = shape,
                             hline = 0, vline = 0,
                             legendPosition = 'right', title = 'PCA bi-plot')

  ## Loading plot
  ans[["loading plot"]] <- plotloadings(p, components=head(p$components, 5),
                                        rangeRetain = 0.01, labSize = 3,
                                        title = 'Loadings plot', axisLabSize = 12,
                                        subtitle = 'PC1, PC2, PC3, PC4, PC5',
                                        caption = 'Top 1% variables',
                                        shape = 24, shapeSizeRange = c(4, 4),
                                        col = c('limegreen', 'black', 'red3'),
                                        legendPosition = 'top',
                                        drawConnectors = FALSE,
                                        returnPlot = TRUE)

  ## eigencorplot
  if(!is.null(metavars)){
    ans[["eigencorplot"]] <- tryCatch(eigencorplot(p, metavars = metavars,
                                                   corFUN = 'pearson',
                                                   corUSE = 'pairwise.complete.obs',
                                                   main = 'PC1-10 clinical correlations'),
                                      error = function(c){NULL})
  }

  return(ans)
}
