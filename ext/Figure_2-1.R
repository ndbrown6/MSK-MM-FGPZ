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

.MMEnv$min_dist = -Inf

registerDoMC(8)

hex_cols = c("#C1272D",
	     "#377EB8",
	     "#01A99D",
	     "#F9ED7D",
	     "#F49C45")

'scientific_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^", scales::scientific_format()(x)), , fixed = TRUE))
}

data("vb=n")
m = data %>% .[["N_Alt"]]
c = data %>% .[["N_Total"]]

# bootstrap symmetric model with actual data vb = n and nb = 1|2|3|4|5|6|7 cell generations
LL = matrix(NA, nrow = 100, ncol = 7)
LL = foreach(i=1:100) %dopar% {
	print(i)
	index = sample(x = 1:length(m), size = rpois(n = 1,lambda = length(m)), replace = TRUE)
	LL = vector(mode = "numeric", length = 7)
	for (j in 1:7) {
		LL[j] = SymmLL(m = m[index], c = c[index], nb = j)$LL
	}
	return(invisible(LL))
}
LL = do.call(rbind, LL)

data_ = dplyr::as_tibble(LL) %>%
	reshape2::melt() %>%
	dplyr::rename(vb = variable,
		      LL = value) %>%
	dplyr::as_tibble() %>%
	dplyr::mutate(vb = gsub("V", "", vb, fixed=TRUE)) %>%
	readr::type_convert() %>%
	dplyr::mutate(n = rep(1:100, times = 7))

plot_ = data_ %>%
	ggplot(aes(x = vb, y = LL/1E3, group = n)) + 
	geom_line(stat = "identity", alpha = .35, color = "#333333", size = .25) +
	theme_classic() +
	xlab(bquote(atop(" ", nu[b]))) +
	ylab(expression("Log-Likelihood ("%.%10^-3~")")) +
	scale_x_continuous(breaks = 1:7,
			   labels = 1:7) +
	scale_y_continuous(breaks = c(0, -1, -2, -3),
			   labels = c("0", "       -1", "-2", "-3"))

pdf(file = "vb_LL.pdf", height = 4, width = 4)
print(plot_)
dev.off()

LL.0 = vector(mode = "numeric", length = 7)
for (i in 1:7) {
	LL.0[i] = SymmLL(m = m, c = c, nb = i)$LL
}

P.0 = vector(mode = "numeric", length = 6)
for (i in 1:6) {
	P.0[i] = 1 - pchisq(2*(LL.0[i+1]-LL.0[i]),1)
}
P.0 = p.adjust(P.0, method = "bonferroni")

# plot density of symmetric model with actual data vb = n and nb = 3|4|5 cell generations
data_ = dplyr::tibble(vaf = m/c)
n = 1E5
dens = list()
for (i in 3:5) {
	LL = SymmLL(m = m, c = c, nb = i)
	vb = .MMEnv$vb[1:i]
	sim_c = sample(x = c, size = n, replace = TRUE)
	sim_v = rbetabinom(n = n, size = sim_c, prob = sample(x = vb, size = n, prob = LL$r_b/sum(LL$r_b), replace = TRUE))/sim_c
	dens[[i-2]] = density(sim_v, from=0, to=0.3, adjust = 1)
}

plot_ = data_ %>%
	ggplot(aes(x = vaf * 100)) +
	geom_histogram(stat = "bin", bins = 15, fill = "#AEA79F", color = "white", alpha = .45) +
	geom_segment(data = dplyr::tibble(x = 1,
					  xend = 30,
				          y = 0,
					  yend = 0),
		  mapping = aes(x = x, y = y, xend = xend, yend = yend),
		  color = "#AEA79F",
		  alpha = .45,
		  inherit.aes = FALSE) +
	geom_vline(xintercept = .MMEnv$vb[1:5]*100, linetype = 3, size = .25) +
	geom_line(data = dplyr::tibble(x = dens[[1]]$x*100,
				       y = dens[[1]]$y),
		  mapping = aes(x = x, y = y),
		  color = hex_cols[1],
		  alpha = .65,
		  size = 1,
		  inherit.aes = FALSE) +
	geom_line(data = dplyr::tibble(x = dens[[2]]$x*100,
				       y = dens[[2]]$y),
		  mapping = aes(x = x, y = y),
		  color = hex_cols[2],
		  alpha = .65,
		  size = 1,
		  inherit.aes = FALSE) +
	geom_line(data = dplyr::tibble(x = dens[[3]]$x*100,
				       y = dens[[3]]$y),
		  mapping = aes(x = x, y = y),
		  color = hex_cols[5],
		  alpha = .65,
		  size = 1,
		  inherit.aes = FALSE) +
	xlab("\nVAF (%)") +
	ylab("Frequency\n\n") +
	scale_x_continuous(limits = c(1,30)) +
	scale_y_continuous(limits = c(0,25)) +
	theme_classic()

pdf(file = "vb_H.pdf", height = 4, width = 4)
print(plot_)
dev.off()

# plot VAF of actual data vb = n and nb = 5 cell generations using symmetric model
LL = SymmLL(m = m, c = c, nb = 5)

data_ = dplyr::tibble(vb = apply(LL$p_bj, 1, which.max)) %>%
	dplyr::mutate(UUID = paste0(data$Gene_Symbol, " ", data$HGVSp_Short),
		      VAF = data$N_Alt/data$N_Total,
		      N_Alt = data$N_Alt,
		      N_Total = data$N_Total,
		      Variant_Classification = data$Variant_Classification) %>%
	dplyr::arrange(vb, desc(VAF)) %>%
	dplyr::mutate(UUID = factor(UUID, levels = unique(UUID), ordered = TRUE)) %>%
	dplyr::mutate(CI95_Lower = binconf(x = N_Alt, n = N_Total, alpha = .05)[,"Lower"]) %>%
	dplyr::mutate(CI95_Upper = binconf(x = N_Alt, n = N_Total, alpha = .05)[,"Upper"])

plot_ = data_ %>%
	dplyr::mutate(Variant_Classification = case_when(
		Variant_Classification == "Frame_Shift_Del" ~ "Frame Shift Deletion",
		Variant_Classification == "Frame_Shift_Ins" ~ "Frame Shift Insertion",
		Variant_Classification == "Missense_Mutation" ~ "Missense Mutation",
		Variant_Classification == "Nonsense_Mutation" ~ "Nonsense Mutation",
		Variant_Classification == "Splice_Site" ~ "Splice Site"
	)) %>%
	ggplot(aes(x = UUID, y = 100*VAF, ymin = CI95_Lower*100, ymax = CI95_Upper*100, fill = factor(vb), shape = Variant_Classification)) +
	geom_hline(yintercept = .MMEnv$vb[1:5]*100, colour = "black", linetype = 3, size = .25) +
	geom_pointrange(show.legend = FALSE) +
	geom_point(stat = "identity", size = 2) +
	scale_fill_manual(values = hex_cols) +
	scale_shape_manual(values = c("Frame Shift Deletion" = 21,
				      "Frame Shift Insertion" = 22,
				      "Missense Mutation" = 23,
				      "Nonsense Mutation" = 24,
				      "Splice Site" = 25)) +
	xlab("\n\n") +
	ylab("\nVAF (%)\n\n") +
	scale_y_continuous(breaks = seq(from = 0, to = 25, by = 5),
			   labels = seq(from = 0, to = 25, by = 5)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
	guides(fill = guide_legend(title = bquote(nu[b]), override.aes = list(shape = 21)),
	       shape = guide_legend(title = "Variant Type"))
	
pdf(file = "VAF_by_Variant.pdf", width = 8, height = 6)
print(plot_)
dev.off()
