#' Boxplot for Olink Data
#'
#' Generate a list of boxplot of NPX vs. grouping variable.
#'
#' @param se a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#' @param features \code{character}(n): the protein names to plot.
#' @param x \code{character}(1): the grouping colname in \code{colData(se)}.
#' @importFrom tibble as_tibble
#' @importMethodsFrom SummarizedExperiment assay
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join filter
#' @importFrom ggplot2 ggplot aes_string geom_boxplot geom_jitter theme ggtitle element_blank element_text
#' @importFrom ggsci scale_fill_npg
#' @importFrom cowplot theme_cowplot
#' @importFrom rlang .data
#' @export
#' @return A list of \code{ggplot} objects
#' @author Ge Tan
#' @examples
#' npxFn <- system.file("extdata",
#'                      c("20200507_Inflammation_NPX_1.xlsx",
#'                        "20200625_Inflammation_NPX_2.xlsx"),
#'                      package = "OlinkR")
#' metaFn <- system.file("extdata", "Inflammation_Metadata.xlsx", package = "OlinkR")
#' se <- as_se(read_npx(npxFn, metaFn))
#' features <- c("IL8", "MCP-3")
#' pList <- olink_boxplot(se, features, x = "condition_Factor")
olink_boxplot <- function(se, features, x = NULL){
  toPlot <- as_tibble(assay(se, "npx"), rownames = "OlinkID") %>%
    pivot_longer(cols = -.data$OlinkID, names_to = "SampleID",
                 values_to = "NPX") %>%
    left_join(as_tibble(colData(se), rownames = "SampleID")) %>%
    left_join(as_tibble(rowData(se), rownames = "OlinkID")) %>%
    filter(.data$Assay %in% features)
  pList <- list()
  feature <- features[1]
  for(feature in features){
    toPlot_one <- filter(toPlot, .data$Assay == feature)
    olinkID <- toPlot_one$OlinkID[1]
    p <- ggplot(toPlot_one, aes_string(x, "NPX", fill = x)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2) +
      scale_fill_npg() +
      ggtitle(str_c(feature, olinkID, sep = " ")) +
      theme_cowplot() +
      theme(axis.text.x=element_blank(),
            plot.title = element_text(hjust = 0.5))
    pList[[feature]] <- p
  }
  return(pList)
}
