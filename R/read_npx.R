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
#' @importFrom OlinkAnalyze read_NPX olink_normalization
#' @importFrom readxl read_excel
#' @importFrom dplyr select right_join bind_rows mutate if_else rename_with across
#' @importFrom tidyselect contains all_of ends_with
#' @importFrom stringr str_c
#' @importFrom magrittr %>%
#' @importFrom stats median model.matrix
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#' @export
#' @return A \code{\link[tibble]{tibble}} in long format with following columns:
#' \describe{
#'   \item{SampleID}{The sample name}
#'   \item{OlinkID}{The internal Olink IDs}
#'   \item{UniProt}{The UniProt IDs}
#'   \item{Assay}{The protein name}
#'   \item{MissingFreq}{The percentage of samples that have intensities lower than LOD}
#'   \item{Panel}{The panel name}
#'   \item{Panel_Version}{The panel version}
#'   \item{PlateID}{The unique plate ID}
#'   \item{QC_Warning}{If the splot passes QC}
#'   \item{LOD}{The limit of detection}
#'   \item{NPX}{NPX value}
#'   \item{Normalization}{The intr-plate normalization method}
#' }
#' @author Ge Tan
#' @examples
#' npxFn <- system.file("extdata",
#'                      c("20200507_Inflammation_NPX_1.xlsx",
#'                        "20200625_Inflammation_NPX_2.xlsx"),
#'                      package = "OlinkR")
#' metaFn <- system.file("extdata", "Inflammation_Metadata.xlsx", package = "OlinkR")
#' npx <- read_npx(npxFn, metaFn)
#'
#' npxFn <- system.file("extdata", "OlinkAnalyze_extdata", "npx_data2.xlsx",
#'                      package = "OlinkR")
#' metaFn <- system.file("extdata", "OlinkAnalyze_extdata" ,
#'                       "npx_data2_meta.xlsx", package = "OlinkR")
#' npx <- read_npx(npxFn, metaFn)
read_npx <- function(npxFn, metaFn, panel = NULL) {
  if(length(npxFn) == 1){
    npx <- read_NPX(npxFn)
  }else{
    # Two cases:
    # 1. different samples from multiple plates
    # 2. some bridging samples from two plates.
    npxList <- map(npxFn, read_NPX)
    ## Sometimes, the NPX files can have different columns from different versions of NPX managers
    overlapCols <- map(npxList, colnames) %>%
      reduce(intersect)
    npxList <- map(npxList, select, all_of(overlapCols))

    repeatedSamples <- map(npxList, function(x){unique(x$SampleID)}) %>%
      reduce(intersect)

    if(length(repeatedSamples) > 0){
      message("Doing reference samples normalisation.")
      if(length(npxList) != 2){
        stop("The reference samples normalisation can only handle two plates.")
        # TODO: if needed, it's possible to extend to more assays.
        #       It's important to choose a reference for all.
      }
      if(length(repeatedSamples) < 8){
        warning("The minimal number of bridging samples for normalisation is 8, you have ", length(repeatedSamples))
      }
      npx <- olink_normalization(df1 = npxList[[1]],
                                 df2 = npxList[[2]],
                                 overlapping_samples_df1 = repeatedSamples)
      npx <- npx[!select(npx, .data$SampleID, .data$OlinkID) %>% duplicated(), ]
    }else{
      message("No enough reference samples. Simply join plates.")
      npx <- bind_rows(npxList)
    }
  }

  npx <- npx %>%
    mutate(Panel = .renamePanels(.data$Panel))

  if (is.null(panel)) {
    panel <- unique(npx$Panel) %>%
      head(1)
  }
  if (length(panel) != 1) {
    stop("This function only reads one panel each time.")
  }
  if (!panel %in% npx$Panel) {
    stop("The panel ", panel, " is not available")
  }
  message("Loading the data from panel: ", panel)

  npx <- npx %>%
    filter(.data$Panel == panel) %>%
    mutate(
      NPX = if_else(is.na(.data$NPX), .data$LOD, .data$NPX),
      isLOD = .data$NPX <= .data$LOD,
      MissingFreq = as.numeric(.data$MissingFreq)
    )
  meta <- read_excel(metaFn)
  meta <- rename_with(meta, make.names)
  meta <- meta %>%
    select(.data$SampleID, ends_with("_Factor"), ends_with("_Numeric")) %>%
    mutate(across(ends_with("_Factor"),
                  function(x){res = make.names(x); res[is.na(x)] = NA; res}))

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

  npx <- npx %>%
    right_join(meta)

  return(npx)
}

#' Coercion to \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#'
#' Coercion to \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param npx \code{\link[tibble]{tibble}} object from \code{\link{read_npx}}.
#' @importFrom magrittr %>%
#' @importFrom dplyr select group_by summarise ungroup distinct
#' @importFrom reshape2 acast
#' @importFrom stringr str_c
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @author Ge Tan
#' @examples
#' npxFn <- system.file("extdata", "OlinkAnalyze_extdata", "npx_data2.xlsx",
#'                      package = "OlinkR")
#' metaFn <- system.file("extdata", "OlinkAnalyze_extdata",
#'                       "npx_data2_meta.xlsx", package = "OlinkR")
#' npx <- read_npx(npxFn, metaFn)
#' se <- as_se(npx)
as_se <- function(npx){
  npxMat <- npx %>%
    select(.data$SampleID, .data$OlinkID, .data$NPX) %>%
    acast(OlinkID ~ SampleID, mean, value.var = "NPX")
  lodMat <- npx %>%
    select(.data$SampleID, .data$OlinkID, .data$LOD) %>%
    acast(OlinkID ~ SampleID, mean, value.var = "LOD")
  colData <- npx %>%
    select(.data$SampleID, .data$PlateID,
           ends_with("_Factor"), ends_with("_Numeric")
           ) %>%
    group_by(.data$SampleID) %>%
    summarise(across(.cols = c("PlateID", ends_with("_Factor")),
                     function(x){str_c(unique(x), collapse = "_")}),
              across(ends_with("_Numeric"), mean, na.rm = TRUE))
  colData <- data.frame(select(colData, -.data$SampleID),
                        row.names = colData$SampleID, check.names = FALSE)
  rowData <- npx %>%
    select(.data$OlinkID, .data$UniProt, .data$Assay,
           .data$Panel, .data$MissingFreq) %>%
    group_by(.data$OlinkID, .data$UniProt, .data$Assay, .data$Panel) %>%
    summarise(MissingFreq = median(.data$MissingFreq)) %>%
    ungroup()
  rowData <- data.frame(select(rowData, -.data$OlinkID),
                        row.names = rowData$OlinkID, check.names = FALSE)

  se <- SummarizedExperiment(
    assays = list(npx = npxMat, npxQCFlag = npxMat > lodMat),
    colData = colData[colnames(npxMat), ],
    rowData = rowData[rownames(npxMat), ]
  )

  return(se)
}

#' List panels from NPX files
#'
#' Read through the given NPX files and return all panels inside
#'
#' @param npxFn Path to NPX files
#' @importFrom OlinkAnalyze read_NPX
#' @importFrom purrr map_dfr
#' @importFrom magrittr %>%
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
  panels <- unique(npx$Panel) %>%
    .renamePanels()
  return(panels)
}

#' Rename panel name
#'
#' Rename panel name to be consistent between different years/batches.
#' @param panels \code{character}(n): the panel names
#' @importFrom dplyr recode
#' @author Ge Tan
.renamePanels <- function(panels){
  recode(panels,
         "Olink INFLAMMATION" = "Olink Target 96 Inflammation")
}
