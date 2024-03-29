#' limma analysis for Olink data
#'
#' Given a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object,
#' perform two groups comparison with \code{limma} package.
#'
#' @param se a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param factorCol \code{character}(1): the metadata column in \code{se} that defines the grouping.
#' @param contrasts \code{character}(1): a character string which can be parsed to expressions, specifying contrasts.
#' @param blocking \code{character}(1): the metadata column in \code{se} which serves as a blocking factor.
#' @importMethodsFrom SummarizedExperiment assay rowData colData
#' @importFrom limma makeContrasts lmFit contrasts.fit topTable eBayes
#' @importFrom stringr str_split str_trim str_replace
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr left_join
#' @importFrom magrittr %>%
#' @importFrom rlang has_name
#' @export
#' @return A \code{tibble} object of results from limma analysis
#' @author Ge Tan
#' @examples
#' npxFn <- system.file("extdata", c(
#'   "20200507_Inflammation_NPX_1.xlsx",
#'   "20200625_Inflammation_NPX_2.xlsx"
#' ),
#' package = "OlinkR"
#' )
#' metaFn <- system.file("extdata", "Inflammation_Metadata.xlsx", package = "OlinkR")
#' se <- as_se(read_npx(npxFn, metaFn))
#' tb <- olink_limma(se,
#'   factorCol = "condition_Factor",
#'   contrasts = "Glucose.10mM.Vehicle - Vehicle.Vehicle",
#'   blocking = "Donor_Factor"
#' )
olink_limma <- function(se, factorCol, contrasts, blocking = NULL) {
  selectedGroups <- str_split(contrasts, pattern = "(-|\\(|\\))")[[1]] %>%
    str_trim()
  if (!has_name(colData(se), factorCol)) {
    stop(factorCol, " column doesn't exist.")
  }
  se <- se[ ,!is.na(colData(se)[[factorCol]])]
  if (!is.null(blocking)) {
    if (!has_name(colData(se), blocking)) {
      stop(blocking, " column doesn't exist.")
    }
    se <- se[ ,!is.na(colData(se)[[blocking]])]
    blocking <- factor(colData(se)[[blocking]])
  }

  eset <- assay(se, "npx")
  Treat <- factor(se[[factorCol]])

  if (is.null(blocking)) {
    design <- model.matrix(~ 0 + Treat)
  } else {
    design <- model.matrix(~ 0 + Treat + blocking)
  }
  colnames(design) <- str_replace(colnames(design), "^Treat", "")
  fit <- lmFit(eset, design)

  cm <- makeContrasts(contrasts = contrasts, levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)

  ans <- topTable(fit2, number = Inf, coef = 1)
  ans <- as_tibble(ans, rownames = "OlinkID")
  ans <- ans %>% left_join(as_tibble(rowData(se), rownames = "OlinkID"))

  selectedCols <- list()
  for (j in 1:length(selectedGroups)) {
    selectedCols[[j]] <- which(Treat == selectedGroups[j])
  }
  useProbe <- logical(nrow(se))
  for (j in 1:length(selectedGroups)) {
    useProbe[rowMeans(assay(se, "npxQCFlag")[, selectedCols[[j]], drop = FALSE]) >= 0.5] <- TRUE
  }
  selectedProbs <- tibble("OlinkID" = rownames(se), isPresent = useProbe)
  selectedExprs <- as_tibble(assay(se, "npx")[, unlist(selectedCols)],
    rownames = "OlinkID"
  )
  ans <- ans %>%
    left_join(selectedProbs) %>%
    left_join(selectedExprs)

  return(ans)
}
