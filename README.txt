###################################
#### Scripts and original data ####
###################################

The initial data consists of the matrices of clonal alterations and epistatic factors in various datasets, including TCGA, paired melanomas and precursors (Shain et al, NEJM 2015), colorectal adenomas and carcinomas (Cross et al, Nat Gen 2018), metastases (Robinson et al, Nature 2017).

## Required general data are collected in the misc folder
CNA_drivers_per_tumor_type.tsv : intOGen CNA drivers
genes_annotation.RData : gene length (in base pairs)
cnSelectL.RData : selection strength for TCGA CNAs, computed from segmented integer copy number data
Mutational_drivers_per_tumor_type.tsv : intOGen mutational drivers

The following scripts have to be run in order to produce all data needed for further processing.

1: tcga/extract_parameters.R
2: models/fit_models_2params.R
3: models/fit_models_3params.R
4: models/select_model.R
5: networks/landscape_networks.R

Epistatic factors are calculated using 10,000 random draws and may therefore vary slightly between replicate simulations. The ones used in the original publications are included in the tcga/epiFactMat_* files. They can however be recomputed during the scripts by commenting/uncommenting code in the last section of the tcga/extract_parameters.R script.

The other scripts can be run in no particular order but may however require to retrieve some external data.

###################################
####       Other datasets      ####
###################################

What needs to be downloaded from third party datasets

## Nevi & Melanomas:
melanoma directory
Need to download supplementary data from the original publication
Mutations: http://www.nejm.org/doi/suppl/10.1056/NEJMoa1502583/suppl_file/nejmoa1502583_appendix_4.xlsx
CNAs: http://www.nejm.org/doi/suppl/10.1056/NEJMoa1502583/suppl_file/nejmoa1502583_appendix_5.xlsx
Script : melanoma_fitness.R
# Change the model if desired (need to comment/uncomment bits of code).

## Coloractal Adenomas & Carcinomas:
adenoma directory


## Metastases (MET500 dataset):
mets directory
The required nature23306-s3.xlsx table for copy number data is available as supplementary table 3 from the original article https://media.nature.com/original/nature-assets/nature/journal/v548/n7667/extref/nature23306-s3.xlsx
The required somatic_v4.csv and cnv_v4.csv tables are available at https://met500.path.med.umich.edu/downloadMet500DataSets
sample_info.xlsx and sequencing_info.xlsx are included in this archive and have been extracted from Supp Tables 1 and 2 (respectively) from the original manuscript's supplementary information in PDF format (https://media.nature.com/original/nature-assets/nature/journal/v548/n7667/extref/nature23306-s1.pdf)


## dN/dS data from Zapata et al:
dnds directory
The required 13059_2018_1434_MOESM2_ESM.txt table of pan-cancer dN/dS per gene can been retrieved from the Additional File 2 from the original article: https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-018-1434-0/MediaObjects/13059_2018_1434_MOESM2_ESM.txt
