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
		  dplyr::left_join(cancer_germ_layer, by = "Case_ID") %>%
		  dplyr::filter(Is_Tissue_Tumor_Normal == "Normal") %>%
		  dplyr::select(Case_ID, Germ_Layer_v2, cell_asymmetry, `VAF_%`) %>%
		  dplyr::mutate(Case_ID = paste0("MOS", Case_ID)) %>%
		  base::as.data.frame()

diverse_index = diversity(data = mosaic_variants %>% dplyr::select(-`VAF_%`), type = c("gini-simpson", "entropy"))
pander::pander(diverse_index[,c("gini.simpson", "entropy")])

data_ = dplyr::as_tibble(mosaic_variants) %>%
	dplyr::left_join(dplyr::tibble(Case_ID = rownames(diverse_index),
				       Gini_Simpson = diverse_index[,"gini.simpson"]), by = "Case_ID")

plot_ = data_ %>%
	ggplot(aes(x = Gini_Simpson, y = cell_asymmetry, fill = Case_ID, shape = Germ_Layer_v2)) +
	geom_point(stat = "identity") +
	scale_shape_manual(values = 21:24) +
	xlab("\nGini-Simpson Index\n") +
	ylab("\nCell division asymmetry\n\n") +
	theme_classic() +
	guides(fill = guide_legend(title = "Patient", override.aes = list(shape = 21), order = 2),
	       shape = guide_legend(title = "Germ layer", order = 1))

pdf(file = "vb_gini.pdf", width = 6, height = 5)
print(plot_)
dev.off()

data_ = dplyr::as_tibble(mosaic_variants) %>%
	dplyr::left_join(dplyr::tibble(Case_ID = rownames(diverse_index),
				       Gini_Simpson = diverse_index[,"entropy"]), by = "Case_ID")

plot_ = data_ %>%
	ggplot(aes(x = Gini_Simpson, y = cell_asymmetry, fill = Case_ID, shape = Germ_Layer_v2)) +
	geom_point(stat = "identity") +
	scale_shape_manual(values = 21:24) +
	xlab("\nShannon Index\n") +
	ylab("\nCell division asymmetry\n\n") +
	theme_classic() +
	guides(fill = guide_legend(title = "Patient", override.aes = list(shape = 21), order = 2),
	       shape = guide_legend(title = "Germ layer", order = 1))

pdf(file = "vb_shannon.pdf", width = 6, height = 5)
print(plot_)
dev.off()

smry_ = mosaic_variants %>%
	dplyr::group_by(Case_ID) %>%	
	dplyr::summarize(sigma = var(cell_asymmetry),
			 mu = mean(cell_asymmetry))

data_ = dplyr::as_tibble(mosaic_variants) %>%
	dplyr::left_join(smry_, by = "Case_ID")

plot_ = data_ %>%
	ggplot(aes(x = sigma, y = cell_asymmetry, fill = Case_ID, shape = Germ_Layer_v2)) +
	geom_point(stat = "identity") +
	scale_shape_manual(values = 21:24) +
	xlab(expression(sigma^2)) +
	ylab("\nCell division asymmetry\n\n") +
	theme_classic() +
	guides(fill = guide_legend(title = "Patient", override.aes = list(shape = 21), order = 2),
	       shape = guide_legend(title = "Germ layer", order = 1))

pdf(file = "vb_sigma.pdf", width = 6, height = 5)
print(plot_)
dev.off()

plot_ = data_ %>%
	ggplot(aes(x = mu, y = cell_asymmetry, fill = Case_ID, shape = Germ_Layer_v2)) +
	geom_point(stat = "identity") +
	scale_shape_manual(values = 21:24) +
	xlab(expression(mu)) +
	ylab("\nCell division asymmetry\n\n") +
	theme_classic() +
	guides(fill = guide_legend(title = "Patient", override.aes = list(shape = 21), order = 2),
	       shape = guide_legend(title = "Germ layer", order = 1))

pdf(file = "vb_mu.pdf", width = 6, height = 5)
print(plot_)
dev.off()

smry_ = mosaic_variants %>%
	dplyr::group_by(Case_ID) %>%	
	dplyr::summarize(sigma = var(`VAF_%`),
			 mu = mean(`VAF_%`))

data_ = dplyr::as_tibble(mosaic_variants) %>%
	dplyr::left_join(smry_, by = "Case_ID")

plot_ = data_ %>%
	ggplot(aes(x = sigma, y = `VAF_%`, fill = Case_ID, shape = Germ_Layer_v2)) +
	geom_point(stat = "identity") +
	scale_shape_manual(values = 21:24) +
	xlab(expression(sigma^2)) +
	ylab("\nVAF (%)\n\n") +
	theme_classic() +
	guides(fill = guide_legend(title = "Patient", override.aes = list(shape = 21), order = 2),
	       shape = guide_legend(title = "Germ layer", order = 1))

pdf(file = "vb_sigma.pdf", width = 6, height = 5)
print(plot_)
dev.off()

plot_ = data_ %>%
	ggplot(aes(x = mu, y = `VAF_%`, fill = Case_ID, shape = Germ_Layer_v2)) +
	geom_point(stat = "identity") +
	scale_shape_manual(values = 21:24) +
	xlab(expression(mu)) +
	ylab("\nVAF (%)\n\n") +
	theme_classic() +
	guides(fill = guide_legend(title = "Patient", override.aes = list(shape = 21), order = 2),
	       shape = guide_legend(title = "Germ layer", order = 1))

pdf(file = "vb_mu.pdf", width = 6, height = 5)
print(plot_)
dev.off()
