---
title: "mouse ortholog"
author: "Ming Tang (edited by Daniel Jeffries for stickleback vs human)"
date: '2024-01-11'
R.version: "4.4.3 (2025-02-28)"
---

#https://bioconductor.org/packages/release/bioc/html/biomaRt.html


```{r}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("XML")


library(dplyr)
library(biomaRt)

listMarts() ## list the databases in ENSEMBL

ensembl <- useEnsembl(biomart = "ensembl") ## connect to the ensemble BioMart

listDatasets(mart = ensembl) ## list all the datasets in the BioMart chosen

##############################################
### HUMAN to STICKLEBACK ###
##############################################

## get human genes with a stickleback ortholog

human<- useMart('ensembl', dataset = "hsapiens_gene_ensembl")

human_attributes<-  c("ensembl_gene_id", "external_gene_name",
                      "gaculeatus_homolog_ensembl_gene", 
                      "gaculeatus_homolog_chromosome",
                      "gaculeatus_homolog_chrom_start",
                      "gaculeatus_homolog_chrom_end",
                      "gaculeatus_homolog_associated_gene_name",
                      "gaculeatus_homolog_orthology_type",
                      "gaculeatus_homolog_perc_id_r1")


listAttributes(human) %>%  ## lists all the attributes of the human dataset relevant for stickleback
  filter(stringr::str_detect(name, "gacu"))

listFilters(human)%>% 
  filter(stringr::str_detect(name, "gaculeatus")) ## check the available filters for gaculeatus

human2stickl_orths <-  getBM(human_attributes, filters="with_gaculeatus_homolog",
                             values=TRUE, mart = human, uniqueRows=T)

length(human2stickl_orths$ensembl_gene_id) 
# >18881 total genes

## How many are one2one orthologs?

human2stickle_one2ones <-  human2stickl_orths %>%
  dplyr::filter(gaculeatus_homolog_orthology_type == "ortholog_one2one")

length(human2stickle_one2ones$ensembl_gene_id)
# >9413

## How many of these are on Chrom 19X ? 

human2stickle_one2ones_chr19 <-  human2stickle_one2ones %>%
  dplyr::filter(gaculeatus_homolog_chromosome == "XIX")

length(human2stickle_one2ones_chr19$ensembl_gene_id)
# >291

### How many on Chrom 12 ? 

human2stickle_one2ones_chr12 <-  human2stickle_one2ones %>%
  dplyr::filter(gaculeatus_homolog_chromosome == "XII")

length(human2stickle_one2ones_chr12$ensembl_gene_id)
# >430 

## There are a few fewer one2ones on Chr19 than we would expect. One reason for this could be that the Y is also present 
## in the assembly, so many of the genes on this linkage group will be one2many (human2stickle). 

## How many one2many are there on Chr19? 

human2stickle_one2many <-  human2stickl_orths %>%
  dplyr::filter(gaculeatus_homolog_orthology_type == "ortholog_one2many")

human2stickle_one2many_chr19 <-  human2stickle_one2many %>%
  dplyr::filter(gaculeatus_homolog_chromosome == "XIX")

length(human2stickle_one2many_chr19$ensembl_gene_id)

# Chr19 - 463

## How does this compare to other chroms? 

human2stickle_one2many_chr_test <-  human2stickle_one2many %>% ## just changing the chrom number one by one
  dplyr::filter(gaculeatus_homolog_chromosome == "XX")

length(human2stickle_one2many_chr_test$ensembl_gene_id)

## random selection of chroms: 

# Chr1 - 397
# Chr4 - 305
# Chr7 - 384
# Chr10 - 270
# Chr13 - 283
# Chr16 - 203
# Chr20 - 223

## Chrom 19Y 
chr19Y_one2ones <-  one2ones %>%
  dplyr::filter(gaculeatus_homolog_chromosome == "Y")

length(chr19Y_one2ones$ensembl_gene_id)



##############################################
### STICKLEBACK to HUMAN ###
##############################################

## get stickleback genes with a human ortholog

stickle <- useMart('ensembl', dataset = "gaculeatus_gene_ensembl")

listAttributes(stickle) %>%  ## search the (large) attribute list for stickleback
  filter(stringr::str_detect(name, "type"))

listAttributes(stickle) %>%  ## search the (large) attribute list for stickleback
  filter(stringr::str_detect(name, "sapien"))

stickle_attributes<-  c("ensembl_gene_id", "external_gene_name", 
                        "chromosome_name", "start_position", "end_position",
                        "gaculeatus_paralog_chromosome", "gaculeatus_paralog_chrom_start",
                        "gaculeatus_paralog_chrom_end", "gaculeatus_paralog_perc_id_r1",
                        "hsapiens_homolog_ensembl_gene", 
                        "hsapiens_homolog_chromosome",
                        "hsapiens_homolog_chrom_start",
                        "hsapiens_homolog_chrom_end",
                        "hsapiens_homolog_associated_gene_name",
                        "hsapiens_homolog_orthology_type",
                        "hsapiens_homolog_perc_id_r1")

## So I cant get a dataset with orthology info AND protein coding genes AND paralogy. I get the error: 
## "Attributes from multiple attribute pages are not allowed"

# So I will make several dataframes, which each contain a set of genes / info that I can then use later to narrow down the final gene set. These will be:

# 1) List of stickleback protein coding genes
# 2) List of stickleback gene paralogs
# 3) List of stickleback > human orthologs. 

## We can then find protein coding genes with paralogs on X and Y that we trust and with orthology. 

## 1) Stickleback protein coding

stickle_protein_coding_attributes<-  c("ensembl_gene_id", "chromosome_name", 
                                       "start_position", "end_position", 
                                       "gene_biotype")




stickle_all_genes <-  getBM(stickle_protein_coding_attributes, 
                                       values=TRUE, 
                                       mart = stickle, 
                                       uniqueRows=T)

stickle_protein_coding_genes <-  stickle_all_genes %>%
  dplyr::filter(gene_biotype == "protein_coding")

length(stickle_protein_coding_genes$ensembl_gene_id)
# > 22341 protein coding genes. 

## 2) Stickleback paralogs

listAttributes(stickle) %>%  ## search the (large) attribute list for stickleback
  filter(stringr::str_detect(name, "paralog"))


stickle_paralog_attributes <-  c("ensembl_gene_id", "chromosome_name", 
                                 "gaculeatus_paralog_ensembl_gene", "gaculeatus_paralog_chromosome",
                                 "gaculeatus_paralog_chrom_start", "gaculeatus_paralog_chrom_end", 
                                 "gaculeatus_paralog_perc_id_r1")


#"start_position", "end_position", 
#"gaculeatus_paralog_chromosome", "gaculeatus_paralog_chrom_start",
#"gaculeatus_paralog_chrom_end",

stickle_paralogs_gene_info <-  getBM(stickle_paralog_attributes, 
                                     values=TRUE, 
                                     mart = stickle, 
                                     uniqueRows=T)
stickle_paralogs_gene_info %>%
  filter(chromosome_name == "XIX", gaculeatus_paralog_perc_id_r1>90)

## Ok reassuringly if I filter for genes on XIX with paralog identity of >90% I get mostly the Y genes! So these are the same genes 


# 3) Stickleback > human orthologs

stickle_human_orthology_attributes<-  c("ensembl_gene_id", "chromosome_name", 
                                        "start_position", "end_position",
                                        "hsapiens_homolog_ensembl_gene", 
                                        "hsapiens_homolog_chromosome",
                                        "hsapiens_homolog_chrom_start",
                                        "hsapiens_homolog_chrom_end",
                                        "hsapiens_homolog_associated_gene_name",
                                        "hsapiens_homolog_orthology_type",
                                        "hsapiens_homolog_perc_id_r1")

stickle2humans_orths <-  getBM(stickle_human_orthology_attributes, 
                               filters="with_hsapiens_homolog",
                               values=TRUE, 
                               mart = stickle,
                               uniqueRows=F)

length(stickle2humans_orths$ensembl_gene_id) 
# >18881 orthologs - pretty good. 

## Output these dataframes
#1)
write.csv(stickle_protein_coding_genes, 
          "~/Data_temp/Stickleback/Stickle_Y_comp/data/orthology/stickleback_UGA_v5_protein_coding_genes.csv")

#2) 
write.csv(stickle_paralogs_gene_info, 
          "~/Data_temp/Stickleback/Stickle_Y_comp/data/orthology/stickleback_UGA_v5_paralogs.csv")

#3) 
write.csv(stickle2humans_orths, 
          "~/Data_temp/Stickleback/Stickle_Y_comp/data/orthology/stickleback_UGA_v5_HG38_orthologs.csv")

####################################################################################
############ Filtering genes for haploinsufficiency analsyses ######################
####################################################################################


## First find one2one orthologs on ChrXIX
stickle2humans_chr19_one2ones <- stickle2humans_orths %>%
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one",
         chromosome_name == "XIX")
names(stickle2humans_chr19_one2ones)

stickle2humans_orths %>% ## note that we don't want these for the HI analysis, but just to show there are about 10 of them 
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one",
         chromosome_name == "Y")

## Next find paralogs on XIX and Y with good identity 

stickle_chr19_XY_paralogs <- stickle_paralogs_gene_info %>%
  filter(gaculeatus_paralog_perc_id_r1>90,
         chromosome_name == "XIX",
         gaculeatus_paralog_chromosome == "Y")
names(stickle_chr19_XY_paralogs)


## Next find one2many orthologs on ChrXIX

stickle2humans_chr19_one2many_all <- stickle2humans_orths %>%
  filter(hsapiens_homolog_orthology_type == "ortholog_one2many",
         chromosome_name == "XIX")

## Now filter the one2manys for those that are good identity paralogs with chrY. I also remove duplicates of the same stickleback gene
## there were only a handful and we don't care about the paralogy, just the orthology. 

stickle2humans_chr19_one2many_filtered <- stickle2humans_chr19_one2many_all %>%
  filter(ensembl_gene_id %in% stickle_chr19_XY_paralogs$ensembl_gene_id,
         ensembl_gene_id %in% stickle_protein_coding_genes$ensembl_gene_id) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

stickle2humans_chr19_one2many_filtered <- stickle_chr19_XY_paralogs %>%
  filter(ensembl_gene_id %in% stickle2humans_chr19_one2many_all$ensembl_gene_id,
         ensembl_gene_id %in% stickle_protein_coding_genes$ensembl_gene_id) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

length(stickle2humans_chr19_one2many_filtered$ensembl_gene_id)

## Now make the Chr19 final gene set (combininensembl_gene_id## Now make the Chr19 final gene set (combining one2ones with the kept one2many)

## this is the basic for all genes (rows)
stickle_gene_info <- stickle_protein_coding_genes %>%
  filter(ensembl_gene_id %in% kept_genes)
names(stickle_gene_info)

## get the info for the paralog genes with human orthology. 

paralog_gene_info <- stickle_protein_coding_genes %>%
  filter(ensembl_gene_id %in% stickle2humans_chr19_one2many_filtered$ensembl_gene_id)

length(paralog_gene_info$ensembl_gene_id)
length(stickle2humans_chr19_one2many_filtered$ensembl_gene_id)

stickle2humans_chr19_one2many_filtered_with_gene_info <- merge(paralog_gene_info,
                                                         stickle2humans_chr19_one2many_filtered)

## I am still missing the orthology info for the one2many genes. 

paralog_ortholog_info <- stickle2humans_chr19_one2many_all %>%
  filter(ensembl_gene_id %in% stickle2humans_chr19_one2many_filtered_with_gene_info$ensembl_gene_id)

length(unique(paralog_ortholog_info$ensembl_gene_id))
length(stickle2humans_chr19_one2many_filtered_with_gene_info$ensembl_gene_id)

## ok now combine these

stickle2humans_chr19_one2many_filtered_with_gene_and_paralog_info <- merge(paralog_ortholog_info,
                                                                           stickle2humans_chr19_one2many_filtered_with_gene_info,
                                                                           all = TRUE)

## So then now I need to combine the one2ones with the kept one2manys with good paralogs

length(stickle2humans_chr19_one2many_filtered_with_gene_and_paralog_info$ensembl_gene_id)
length(stickle2humans_chr19_one2ones$ensembl_gene_id)

names(stickle2humans_chr19_one2many_filtered_with_gene_and_paralog_info)
names(stickle2humans_chr19_one2ones)


stickle2_humans_all_orths <- merge(stickle2humans_chr19_one2ones,
                                   stickle2humans_chr19_one2many_filtered_with_gene_and_paralog_info,
                                   all = TRUE)

stickle2_humans_all_orths$gene_biotype = "protein_coding" ## this got missed from the "one2ones" but not worth specifically filtering and merging just for this - all genes are protein coding so just adding manually here
stickle2_humans_all_orths

length(stickle2_humans_all_orths$ensembl_gene_id)

write.csv(stickle2_humans_all_orths, 
          "~/Data_temp/Stickleback/Stickle_Y_comp/results/BiomaRt/stickle_2_humans_all_orths.csv",
          row.names = F)

