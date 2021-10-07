####################################################
##### Installing and loading required packages
####################################################



if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("gam")) {
  install.packages("glmnet", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("rdist")) {
  install.packages("doMC", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("dynamicTreeCut")) {
  install.packages("StatMatch", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("circlize")) {
  install.packages("Rtsne", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("rlist")) {
  install.packages("fpc", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("gridExtra")) {
  install.packages("GA", dependencies = TRUE, repos="http://cran.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}

if (!require("Rsamtools")) {BiocManager::install("Rsamtools")}

if (!require("Biostrings")) {BiocManager::install("Biostrings")}



if (!require("GenomicRanges")) {BiocManager::install("GenomicRanges")}

if (!require("regioneR")) {BiocManager::install("regionR")}

if (!require("IRanges")) {BiocManager::install("IRanges")}


if (!require("karyoploteR")) {BiocManager::install("regioneR")}



#Load libraries
library(ggplot2)
library(data.table)
library(gam)
library(rdist)
library(dynamicTreeCut)
library(circlize)
library(rlist)
library(gridExtra)
library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(IRanges)
library(karyoploteR)

