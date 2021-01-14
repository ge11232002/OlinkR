#' Produce files for webgestalt website
#'
#' Given a \code{tibble} object of the results returned from \code{olink_limma} function, generate the files for upload to webgestalt website.
#'
#' @param x A \code{tibble} object from \code{olink_limma} function.
#' @param pvalue The p-value to filter for significant proteins. The default is 0.01.
#' @param log2FC The log2 fold change cut-off to filter for significant proteins. The default is 0, i.e. no filtering.
#' @param dir The output folder. The default is current working directory.
#' @importFrom readr write_lines write_tsv
#' @importFrom dplyr filter pull select
#' @export
#' @return An invisible character object with exported file names.
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
#' files <- webgestalt_prep(se)
#' file.remove(files)
webgestalt_prep <- function(x, pvalue=0.01, log2FC=0, dir="."){
  dir.create(dir, recursive = TRUE)

# ORA files ---------------------------------------------------------------
  tb_ora <- x %>% filter(logFC > log2FC, P.Value <= pvalue) %>%
    pull(UniProt)
  write_lines(tb_ora, file = file.path(dir, "webgestalt_upregulated_ora.txt"))

  tb_ora <- x %>% filter(logFC < -log2FC, P.Value <= pvalue) %>%
    pull(UniProt)
  write_lines(tb_ora, file = file.path(dir, "webgestalt_downregulated_ora.txt"))

# GSEA --------------------------------------------------------------------
  tb_gsea <- x %>% select(UniProt, logFC)
  write_tsv(tb_gsea, file = file.path(dir, "webgestalt_gsea.rnk"),
            col_names = FALSE)

# Background --------------------------------------------------------------
  tb_background <- x %>% pull(UniProt)
  write_lines(tb_background, file = file.path(dir, "webgestalt_background.txt"))

  invisible(file.path(dir, c("webgestalt_upregulated_ora.txt",
                             "webgestalt_downregulated_ora.txt",
                             "webgestalt_gsea.rnk", "webgestalt_background.txt")))
}
