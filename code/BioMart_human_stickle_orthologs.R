install.packages("BiocManager")
BiocManager::install("biomaRt")


---
  title: "mouse ortholog"
author: "Ming Tang (edited by Daniel Jeffries for stickleback vs human)"
date: '2024-01-11'
R.version: "4.4.3 (2025-02-28)"
output: html_document
editor_options: 
  chunk_output_type: console
---
  
  Get mouse orthologs for human

https://bioconductor.org/packages/release/bioc/html/biomaRt.html



```{r}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("XML")


library(dplyr)
library(biomaRt)

listMarts() ## list the databases in ENSEMBL

ensemble <- useEnsembl(biomart = "ensembl") ## connect to the ensemble BioMart

listDatasets(mart = ensemble) ## list all the datasets in the BioMart chosen

## get human genes with a stickleback ortholog

human<- useMart('ensembl', dataset = "hsapiens_gene_ensembl")

attributes<-  c("ensembl_gene_id", "external_gene_name",
                "gaculeatus_homolog_ensembl_gene", 
                "gaculeatus_homolog_chromosome",
                "gaculeatus_homolog_chrom_start",
                "gaculeatus_homolog_chrom_end",
                "gaculeatus_homolog_associated_gene_name",
                "gaculeatus_homolog_orthology_type",
                "gaculeatus_homolog_perc_id_r1")

listAttributes(human) %>%  ## lists all the attributes of the dataset (there are a lot)
  head()

listAttributes(human) %>%  ## search the (large) attribute list for stickleback
  filter(stringr::str_detect(name, "gacu"))

listAttributes(human) %>%  ## search the (large) attribute list for stickleback
  filter(stringr::str_detect(name, "sapien"))

listAttributes(human)  %>% head()



orth.stickle_nonUniq <-  getBM(attributes, filters="with_gaculeatus_homolog",
                    values=TRUE, mart = human, uniqueRows=T)

listFilters(human)%>% 
  filter(stringr::str_detect(name, "gaculeatus"))

length(orth.stickle$ensembl_gene_id)

```
## get stickleback genes with a human ortholog

#stickle <- useMart('ensembl', dataset = "gaculeatus_gene_ensembl")

#attributes<-  c("ensembl_gene_id", "external_gene_name",
#                "gaculeatus_homolog_ensembl_gene", 
#                "gaculeatus_homolog_chromosome",
#                "gaculeatus_homolog_chrom_start",
#                "gaculeatus_homolog_chrom_end",
#                "gaculeatus_homolog_associated_gene_name",
#                "gaculeatus_homolog_orthology_type",
#                "gaculeatus_homolog_perc_id_r1")

#orth.stickle_nonUniq <-  getBM(attributes, filters="with_gaculeatus_homolog",
#                               values=TRUE, mart = human, uniqueRows=F)

#unique(orth.stickle$gaculeatus_homolog_orthology_type)


#listFilters(human)%>% head()
#listFilters(human)%>% 
#  filter(stringr::str_detect(name, "gaculeatus"))

#length(orth.stickle$ensembl_gene_id)



## Filter the one2one orthologs

```{r}
human_to_stickle_one2ones <-  orth.stickle %>%
  dplyr::filter(gaculeatus_homolog_orthology_type == "ortholog_one2one")

length(one2ones$ensembl_gene_id)

one2ones_nonUnique <-  orth.stickle_nonUniq %>%
  dplyr::filter(gaculeatus_homolog_orthology_type == "ortholog_one2one")

length(one2ones_nonUnique$ensembl_gene_id)

### Chrom 19 

chr19_one2ones <-  one2ones %>%
  dplyr::filter(gaculeatus_homolog_chromosome == "XIX")

length(chr19_one2ones$ensembl_gene_id)


chr19_one2ones <-  one2ones %>%
  dplyr::filter(gaculeatus_homolog_chromosome == "XIX")

length(chr19_one2ones$ensembl_gene_id)

### Chrom 12 

chr12_one2ones <-  one2ones %>%
  dplyr::filter(gaculeatus_homolog_chromosome == "XII")

length(chr12_one2ones$ensembl_gene_id)

## Chrom 19Y 
chr19Y_one2ones <-  one2ones %>%
  dplyr::filter(gaculeatus_homolog_chromosome == "Y")

length(chr19Y_one2ones$ensembl_gene_id)

## one2many
one2many <-  orth.stickle %>%
  dplyr::filter(gaculeatus_homolog_orthology_type == "ortholog_one2many")

length(one2many$ensembl_gene_id) ## expect ~300

## Chr19

Chr19_one2many <-  one2many %>%
  dplyr::filter(gaculeatus_homolog_chromosome == "XIX")

length(Chr19_one2many$ensembl_gene_id)

## Chr12

Chr12_one2many <-  one2many %>%
  dplyr::filter(gaculeatus_homolog_chromosome == "XII")

length(Chr12_one2many$ensembl_gene_id)


## Chr19Y 
Chr19Y_one2many <-  one2many %>%
  dplyr::filter(gaculeatus_homolog_chromosome == "Y")

length(Chr19Y_one2many)

```