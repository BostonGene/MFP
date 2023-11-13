install.packages("BiocManager")

path_to_download_packages<-.libPaths()

BiocManager::install(c("affy", "annotate", "gcrma"), dependencies=TRUE, lib=path_to_download_packages[1])