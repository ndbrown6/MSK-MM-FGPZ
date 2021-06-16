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
		  dplyr::filter(Is_Mosaic_Yes_No == "Yes") %>%
		  dplyr::mutate(cell_asymmetry = `VAF_%`/(vb*100))

cancer_germ_layer = mosaic_variants %>%
		    dplyr::filter(Is_Mosaic_Yes_No == "Yes",
				  Is_Tissue_Tumor_Normal == "Tumor") %>%
		    dplyr::select(Case_ID, Cancer_Germ_Layer_v1 = Germ_Layer_v1, Cancer_Germ_Layer_v2 = Germ_Layer_v2) %>%
		    dplyr::filter(!duplicated(paste0(Case_ID, ":", Cancer_Germ_Layer_v2)))

mosaic_variants = mosaic_variants %>%
		  dplyr::left_join(cancer_germ_layer, by = "Case_ID")

data_ = mosaic_variants %>%
	dplyr::mutate(y = cell_asymmetry) %>%
	dplyr::mutate(x = case_when(
		Germ_Layer_v2 == Cancer_Germ_Layer_v2 ~ Is_Tissue_Tumor_Normal,
		TRUE ~ "Other"
	)) %>%
	dplyr::mutate(x = case_when(
		x == "Tumor" ~ "Tumor",
		x == "Normal" ~ "Tumor-matched\ngerm layer",
		x == "Other" ~ "Other non-tumor\ngerm layers"
	)) %>%
	dplyr::mutate(x = factor(x, levels = c("Other non-tumor\ngerm layers", "Tumor-matched\ngerm layer", "Tumor"), ordered = TRUE))

plot_ = data_ %>%
	ggplot(aes(x = x, y = y)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA) +
	geom_jitter(data = data_, mapping = aes(x = x, y = y, fill = factor(Case_ID), shape = factor(Germ_Layer_v2)),
		    stat = "identity", width = .10, height = 0, size = 2.5, alpha = .85, inherit.aes = FALSE) +
	scale_shape_manual(values = 21:24) +
	xlab("") +
	ylab("\nCell division asymmetry\n\n") +
	scale_y_sqrt(breaks = c(1, 5, 10, 15, 20),
		     labels = c(1, 5, 10, 15, 20)) +
	theme_classic() +
	guides(fill = guide_legend(title = "Patient", override.aes = list(shape = 21)),
	       shape = guide_legend(title = "Germ layer")) +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE),
		    y_position = sqrt(10),
		    comparison = list(c("Other non-tumor\ngerm layers", "Tumor-matched\ngerm layer")))

pdf(file = "vb_all.pdf", width = 6, height = 5)
print(plot_)
dev.off()
