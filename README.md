# OlinkR
Data analysis package for Olink data

## Installation
`OlinkR` package has an non-CRAN/Bioconductor dependency `OlinkAnalyze` on github,
and several other Bioconductor packages.
To install them,
```R
remotes::install_github(repo ='Olink-Proteomics/OlinkRPackage/OlinkAnalyze',
                        ref = "main", build_vignettes = FALSE)
BiocManager::install(c("PCAtools", "scater"))
```

Then install `OlinkR` package.
```R
remotes::install_github("ge11232002/OlinkR", dependencies = TRUE)
```

