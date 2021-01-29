#' Read NPX and metadata file together
#'
#' Read one or more NPX files from Olink NPX manager and
#' one metadata file in excel format.
#'
#' The metadata excel file shall have column \dQuote{SampleID} and other metadata column that ends with \dQuote{_Factor} (mandatory )or \dQuote{_Numeric} (optional).
#' If some samples from NPX don't exist in metadata file, it gives a warning.
#' Only the samples in metadata file are kept.
#' If some samples from metadata file don't exists in NPX file, it gives an error.
#' @param npxFn Path to NPX files
#' @param metaFn Path to metadata excel file
#' @param panel \code{character}(1): the panel to load.
#'              By default, the first panel in the data will be loaded.
#' @importFrom OlinkAnalyze read_NPX
#' @importFrom readxl read_excel
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom reshape2 acast
#' @importFrom dplyr select right_join group_by summarise ungroup bind_rows mutate if_else rename_with across
#' @importFrom tidyselect contains
#' @importFrom stringr str_c
#' @importFrom magrittr %>%
#' @importFrom stats median model.matrix
#' @importFrom purrr map_dfr
#' @export
#' @return A list with two objects: a \code{\link[tibble]{tibble}} in long format and a
#'         \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @author Ge Tan
#' @examples
#' npxFn <- system.file("extdata", c(
#'   "20200507_Inflammation_NPX_1.xlsx",
#'   "20200625_Inflammation_NPX_2.xlsx"
#' ),
#' package = "OlinkR"
#' )
#' metaFn <- system.file("extdata", "Inflammation_Metadata.xlsx", package = "OlinkR")
#' read_npx(npxFn, metaFn)
read_npx <- function(npxFn, metaFn, panel = NULL) {
  npx <- map_dfr(npxFn, read_NPX)
  if (is.null(panel)) {
    panel <- unique(npx$Panel) %>% head(1)
  }
  if (length(panel) != 1) {
    stop("This function only reads one panel each time.")
  }
  if (!panel %in% npx$Panel) {
    stop("The panel ", panel, " is not available")
  }
  message("Loading the data from panel: ", panel)

  npx <- npx %>%
    filter(Panel == panel) %>%
    mutate(
      NPX = if_else(is.na(NPX), LOD, NPX), isLOD = NPX <= LOD,
      MissingFreq = as.numeric(MissingFreq)
    )
  meta <- read_excel(metaFn)
  meta <- rename_with(meta, make.names)
  meta <- meta %>% select(
    "SampleID", ends_with("_Factor"), ends_with("_Numeric")
  ) %>%
    mutate(across(ends_with("_Factor"), function(x){res=make.names(x);res[is.na(x)]=NA;res}))

  if (length(setdiff(npx$SampleID, meta$SampleID)) != 0L) {
    warning(
      "SampleID (",
      str_c(setdiff(npx$SampleID, meta$SampleID), collapse = ", "),
      ") is/are not in metadata."
    )
  }
  if (length(setdiff(meta$SampleID, npx$SampleID)) != 0L) {
    stop(
      "SampleID (",
      str_c(setdiff(meta$SampleID, npx$SampleID), collapse = ", "),
      ") is/are not in NPX data."
    )
  }

  npx <- npx %>% right_join(meta)

  ## Build SummarizedExperiment object ------------------------------
  npxMat <- npx %>%
    select(SampleID, OlinkID, NPX) %>%
    acast(OlinkID ~ SampleID, value.var = "NPX")
  lodMat <- npx %>%
    select(SampleID, OlinkID, LOD) %>%
    acast(OlinkID ~ SampleID, value.var = "LOD")
  colData <- npx %>%
    select(PlateID, colnames(meta)) %>%
    unique()
  colData <- data.frame(
    select(colData, -SampleID),
    row.names = colData$SampleID, check.names = FALSE
  )
  rowData <- npx %>%
    select(OlinkID, UniProt, Assay, Panel, MissingFreq) %>%
    group_by(OlinkID, UniProt, Assay, Panel) %>%
    summarise(MissingFreq = median(MissingFreq)) %>%
    ungroup()
  rowData <- data.frame(
    select(rowData, -OlinkID),
    row.names = rowData$OlinkID, check.names = FALSE
  )

  se <- SummarizedExperiment(
    assays = list(npx = npxMat, npxQCFlag = npxMat > lodMat),
    colData = colData[colnames(npxMat), ],
    rowData = rowData[rownames(npxMat), ]
  )
  return(list(tibble = npx, SummarizedExperiment = se))
}

#' List panels from NPX files
#'
#' Read through the given NPX files and return all panels inside
#'
#' @param npxFn Path to NPX files
#' @importFrom OlinkAnalyze read_NPX
#' @importFrom purrr map_dfr
#' @export
#' @return A \code{character}(n) object.
#' @author Ge Tan
#' @examples
#' npxFn <- system.file("extdata", c(
#'   "2019-01-02_20181145_Akdis_Nasal_NPX_LOD.xlsx"
#' ),
#' package = "OlinkR"
#' )
#' list_panels(npxFn)
list_panels <- function(npxFn) {
  npx <- map_dfr(npxFn, read_NPX)
  panels <- unique(npx$Panel)
  return(panels)
}
