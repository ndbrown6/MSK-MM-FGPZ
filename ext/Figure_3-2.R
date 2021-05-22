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

registerDoMC(8)

hex_cols = c("#C1272D",
	     "#377EB8",
	     "#01A99D",
	     "#F9ED7D",
	     "#F49C45")

'scientific_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^", scales::scientific_format()(x)), , fixed = TRUE))
}

# predict cell generation from actual vb=n data
data("vb=n")
m = data %>% .[["N_Alt"]]
c = data %>% .[["N_Total"]]

LL = MixtureBB(m = m, c = c, nb = 4)
.MMEnv$vb = LL$vb
.MMEnv$max_iter = 1000
LL0 = SymmBB(m = m, c = c, nb = 4)

blood_variants = dplyr::tibble(nb = apply(LL0$p_bj, 1, which.max)) %>%
		 dplyr::mutate(vb = LL$vb[nb]) %>%
		 dplyr::mutate(UUID = paste0(data$Gene_Symbol, " ", data$HGVSp_Short))

# predict cell generation from actual vb=n data
data("vb=all")

all_variants = data %>%
	       dplyr::mutate(`VAF_%` = gsub(pattern = "%", replacement = "", x = `VAF_%`, fixed = TRUE)) %>%
	       readr::type_convert() %>%
	       dplyr::mutate(UUID = paste0(data$Gene_Symbol, " ", data$HGVSp_Short)) %>%
	       dplyr::left_join(blood_variants, by = "UUID") %>%
	       dplyr::filter(Is_Tissue_Tumor_Normal == "Normal",
			     Is_Mosaic_Yes_No == "Yes") %>%
	       dplyr::mutate(cell_asymmetry = `VAF_%`/(vb*100))


plot_ = all_variants %>%
	ggplot(aes(x = Germ_Layer_v2, y = cell_asymmetry)) +
	geom_boxplot(stat = "boxplot") +
	facet_wrap(~Case_ID, scales = "free_x")