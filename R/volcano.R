#' Volcano plot of Olink results
#'
#' Given the \code{tibble} object from \code{olink_limma}, plot a volcano plot.
#' The log2FC is on the x-axis and -log10(P.value) is on y-axis.
#' By default, significant proteins are labeled on the volcano plot.
#'
#' @param tb A \code{tibble} object from \code{olink_limma} function.
#' @param p.value \code{numeric}(1): The p.value cutoff for significance level.
#' @param log2FC \code{numeric}(1): The log2FC cutoff for significance level.
#' @param olinkIds \code{character}(n): OlinkIDs to label on the volcano plot.
#' @importFrom cowplot theme_cowplot
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr mutate case_when filter arrange pull
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_x_continuous geom_hline geom_vline xlab ylab
#' @importFrom rlang .data
#' @export
#' @return A \code{ggplot} object.
#' @author Ge Tan
#' @examples
#' npxFn <- system.file("extdata",
#'                      c("20200507_Inflammation_NPX_1.xlsx",
#'                        "20200625_Inflammation_NPX_2.xlsx"),
#'                      package = "OlinkR")
#' metaFn <- system.file("extdata", "Inflammation_Metadata.xlsx",
#'                       package = "OlinkR")
#' se <- as_se(read_npx(npxFn, metaFn))
#' tb <- olink_limma(se,
#'   factorCol = "condition_Factor",
#'   contrasts = "Glucose.10mM.Vehicle - Vehicle.Vehicle",
#'   blocking = "Donor_Factor"
#' )
#' olink_volcano(tb)
olink_volcano <- function(tb, p.value = 0.05, log2FC = 0, olinkIds = NULL) {
  tb <- tb %>%
    mutate(type = case_when(
      .data$P.Value <= p.value & abs(.data$logFC) >= log2FC & .data$isPresent ~ "Significant",
      (.data$P.Value > p.value | abs(.data$logFC) < log2FC) & .data$isPresent ~ "Non-significant",
      TRUE ~ "Absent")
      )
  if (is.null(olinkIds)) {
    olinkIds <- tb %>%
      filter(.data$type == "Significant") %>%
      pull(.data$OlinkID)
  }

  p <- ggplot(tb %>% arrange(.data$type),
              aes(.data$logFC, -log10(.data$P.Value))) +
    geom_point(aes(col = .data$type)) +
    geom_label_repel(
      data = filter(tb, .data$OlinkID %in% olinkIds),
      aes(label = .data$Assay),
      box.padding = 1,
      segment.color = "grey50", max.overlaps = 40) +
    scale_color_manual(values = c(
      "Significant" = "red", "Non-significant" = "black",
      "Absent" = "gray")) +
    geom_hline(yintercept = -log10(p.value), linetype = "dotted") +
    geom_vline(xintercept = log2FC, linetype = "dotted") +
    geom_vline(xintercept = -log2FC, linetype = "dotted") +
    xlab(expression(paste(log[2], "FC"))) +
    ylab(expression(paste(-log[10], "pvalue"))) +
    scale_x_continuous(limits = c(-1.2 * max(abs(tb$logFC)), 1.2 * max(abs(tb$logFC)))) +
    theme_cowplot()

  return(p)
}
