#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
library('mosaicm')
library('ggplot2')
library('ggridges')
library('parallel')
library('foreach')
library('doMC')
library('Hmisc')
library('copynumber')
library('ggsignif')
library('diverse')

data("vb=all")

mosaic_variants = data %>%
		  dplyr::mutate(Sample_ID = gsub(pattern = "-", replacement = "", x = Sample_ID, fixed = TRUE)) %>%
		  dplyr::mutate(`VAF_%` = gsub(pattern = "%", replacement = "", x = `VAF_%`, fixed = TRUE)) %>%
		  readr::type_convert() %>%
		  dplyr::mutate(UID = paste0(Case_ID, ":", Gene_Symbol, ":", HGVSp_Short)) %>%
		  dplyr::mutate(UUID = paste0(Case_ID, ":", Sample_ID, ":", Gene_Symbol, ":", HGVSp_Short)) %>%
		  dplyr::filter(Is_Mosaic_Yes_No == "Yes") %>%
		  dplyr::filter(Is_Tissue_Tumor_Normal == "Normal") %>%
		  dplyr::select(Case_ID, Germ_Layer_v2, `VAF_%`) %>%
		  dplyr::mutate(Case_ID = paste0("MOS", Case_ID)) %>%
		  base::as.data.frame()

d = diversity(data = mosaic_variants, type = c("gini-simpson", "entropy"))
pander::pander(d[,c("gini.simpson", "entropy")])
