
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXmixomics")

library(LXmixomics)

??LXmixomics
#---------------------
data(gene_data)
data(protein_data)
data(meta_data)

#--------------------

rm(list=ls())


gene_data= "gene_example.xlsx"
protein_data= "protein_example.xlsx "
meta_data="meta_example.xlsx "

#devtools::load_all()

LXmixomics(gene_data,protein_data,meta_data)



