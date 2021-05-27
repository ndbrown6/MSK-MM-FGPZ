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
		 dplyr::mutate(UID = paste0(data$Case_ID, ":", data$Gene_Symbol, ":", data$HGVSp_Short))

# predict cell generation from actual vb=n data
data("vb=all")

mosaic_variants = data %>%
		  dplyr::mutate(Sample_ID = gsub(pattern = "-", replacement = "", x = Sample_ID, fixed = TRUE)) %>%
		  dplyr::mutate(`VAF_%` = gsub(pattern = "%", replacement = "", x = `VAF_%`, fixed = TRUE)) %>%
		  readr::type_convert() %>%
		  dplyr::mutate(UID = paste0(Case_ID, ":", Gene_Symbol, ":", HGVSp_Short)) %>%
		  dplyr::mutate(UUID = paste0(Case_ID, ":", Sample_ID, ":", Gene_Symbol, ":", HGVSp_Short)) %>%
		  dplyr::left_join(blood_variants, by = "UID") %>%
		  dplyr::filter(Is_Mosaic_Yes_No == "Yes")

ampliseq_variants = readr::read_tsv(file = "~/Desktop/MOSAICISM - Mutation file - with Ampliseq- 04-25-21 - with zeros.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		    readr::type_convert() %>%
		    dplyr::mutate(UUID = paste0(Case_ID, ":", Sample_ID, ":", Gene_Symbol, ":", HGVSp_Short)) %>%
		    dplyr::filter(!is.na(AMPLISEQ_MAF)) %>%
		    dplyr::filter(!(Case_ID %in% c(1, 14))) %>%
		    dplyr::select(AMPLISEQ_MAF, AMPLISEQ_ALT_Count, AMPLISEQ_Total_Count, UUID)

mosaic_variants = mosaic_variants %>%
#		  dplyr::left_join(ampliseq_variants, by = "UUID") %>%
#		  dplyr::mutate(AMPLISEQ_MAF = AMPLISEQ_MAF * 100) %>%
#		  dplyr::mutate(`VAF_%` = case_when(
#		  	is.na(AMPLISEQ_MAF) ~ `VAF_%`,
#			TRUE ~ AMPLISEQ_MAF
#		  )) %>%
		  dplyr::mutate(cell_asymmetry = `VAF_%`/(vb*100))


cancer_germ_layer = mosaic_variants %>%
		    dplyr::filter(Is_Mosaic_Yes_No == "Yes",
				  Is_Tissue_Tumor_Normal == "Tumor") %>%
		    dplyr::select(Case_ID, Cancer_Germ_Layer_v1 = Germ_Layer_v1, Cancer_Germ_Layer_v2 = Germ_Layer_v2)


mosaic_variants = mosaic_variants %>%
		  dplyr::left_join(cancer_germ_layer, by = "Case_ID")

plot_ = mosaic_variants %>%
	dplyr::mutate(y = cell_asymmetry) %>%
	dplyr::mutate(x = case_when(
		Is_Tissue_Tumor_Normal == "Tumor" ~ paste0(Is_Tissue_Tumor_Normal, "\n", Germ_Layer_v1),
		Germ_Layer_v1 == Cancer_Germ_Layer_v1 ~ paste0(Is_Tissue_Tumor_Normal, "\n", Germ_Layer_v1),
		TRUE ~ paste0(Is_Tissue_Tumor_Normal, "\n Other germlayers")
	)) %>%
	ggplot(aes(x = x, y = y)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA) +
	geom_jitter(stat = "identity", width = .10, height = 0, shape = 21, size = 2.5, alpha = .75) +
	xlab("\n\nGermlayer\n") +
	ylab("\nCell division asymmetry\n\n") +
	scale_y_sqrt() +
	facet_wrap(~Case_ID, scales = "free")

pdf(file = "vb_all.pdf", width = 14, height = 9)
print(plot_)
dev.off()

