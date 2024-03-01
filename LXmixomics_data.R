
library(openxlsx)

gene_data <- read.xlsx("gene_example.xlsx")
protein_data <- read.xlsx("protein_example.xlsx")
meta_data <- read.xlsx("meta_example.xlsx")

usethis::use_data(gene_data,overwrite = T)
usethis::use_data(protein_data,overwrite = T)
usethis::use_data(meta_data,overwrite = T)

rm(list=ls())

data(gene_data)
data(protein_data)
data(meta_data)
