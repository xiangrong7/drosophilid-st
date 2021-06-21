install.packages(c("dplyr","Seurat","SeuratObject","geomorph",
                   "plotly","RColorBrewer","readobj",
                   "rgl","dashCoreComponents","dash","dashHtmlComponents"))
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

