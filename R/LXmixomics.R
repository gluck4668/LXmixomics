

LXmixomics <- function(gene_data,protein_data,meta_data){

# https://www.jianshu.com/p/734d6fc3ecd6

# Integrative analysis of metabolomics and proteomics reveals amino acid metabolism disorder in sepsis

#---------安装相关R包---------------------------------------------------------
packges_isntall <- function(){

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  installed_packs <- installed.packages()[,1] #已安装的R包
  com_packs <- c("igraph","rgl","ellipse","corpcor","RColorBrewer","plyr",
                 "parallel","dplyr","tidyr","reshape2","methods","matrixStats",
                 "rARPACK","gridExtra","MASS","lattice","ggplot2","openxlsx",
                 "limma","stringr") # 常规R包
  bio_packs <- c("mixOmics") # BiocManager包

  # 未安装的R包
  not_com <- com_packs[!com_packs %in% installed_packs]
  not_bio <- bio_packs[!bio_packs %in% installed_packs]

  if(length(not_com)>0){
    instll_com <- function(i){install.packages(i,ask=F,update=F)}
    sapply(not_packs,instll_com,simplify = T)
  }

  if(length(not_bio)>0){
    instll_bio <- function(i){BiocManager::install(i,ask=F,update=F)}
    sapply(not_bio,instll_bio,simplify = T)
  }


  packs <- c(com_packs,bio_packs)
  lib_fun <- function(i){library(i,character.only = T)}
  sapply(packs,lib_fun,simplify = T)

}

packges_isntall()

#-----------------------------------------------------------------------------

gene_data=trimws(gene_data)
protein_data=trimws(protein_data)
meta_data=trimws(meta_data)

if(length(gene_data)==0)
  gene_data=""
if(length(protein_data)==0)
  protein_data=""
if(length(meta_data)==0)
  meta_data=""

if(str_length(gene_data)<1)
   LXmixomics_pro_meta(protein_data,meta_data) else
     {if(str_length(protein_data)<1)
       LXmixomics_gene_meta(gene_data,meta_data) else
         if(str_length(meta_data)<1)
           LXmixomics_gene_pro(gene_data,protein_data) else
             LXmixomics_3omics(gene_data,protein_data,meta_data)
     }


}
