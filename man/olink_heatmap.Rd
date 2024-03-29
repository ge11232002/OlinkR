% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap.R
\name{olink_heatmap}
\alias{olink_heatmap}
\title{Heatmap of Olink results}
\usage{
olink_heatmap(tb, se, p.value = 0.05, log2FC = 0, ...)
}
\arguments{
\item{tb}{A \code{tibble} object from \code{olink_limma} function.}

\item{se}{A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.}

\item{p.value}{\code{numeric}(1): The p.value cutoff for significance level.}

\item{log2FC}{\code{numeric}(1): The log2FC cutoff for significance level.}

\item{...}{graphical parameter passed to \code{pheatmap}.}
}
\value{
An invisible \code{pheatmap} object.
}
\description{
Given the \code{tibble} object from \code{olink_limma}, plot a heatmap of significant
proteins.
}
\examples{
npxFn <- system.file("extdata",
                     c("20200507_Inflammation_NPX_1.xlsx",
                       "20200625_Inflammation_NPX_2.xlsx"),
                     package = "OlinkR")
metaFn <- system.file("extdata", "Inflammation_Metadata.xlsx", package = "OlinkR")
se <- as_se(read_npx(npxFn, metaFn))
tb <- olink_limma(se,
                  factorCol = "condition_Factor",
                  contrasts = "Glucose.10mM.Vehicle - Vehicle.Vehicle",
                  blocking = "Donor_Factor")
olink_heatmap(tb, se)
}
\author{
Ge Tan
}
