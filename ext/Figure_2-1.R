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

registerDoMC(12)

hex_cols = c("#e41a1c",
	     "#377eb8",
	     "#4daf4a",
	     "#984ea3",
	     "#ff7f00")

'scientific_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^", scales::scientific_format()(x)), , fixed = TRUE))
}

data("vb=n")
m = data %>% .[["N_Alt"]]
c = data %>% .[["N_Total"]]

# test symmetric model with actual data vb = n and nb = 1|2|3|4|5|6|7 cell generations
LL = matrix(NA, nrow = 100, ncol = 7)
LL = foreach(i=1:100) %dopar% {
	index = sample(x = 1:length(m), size = length(m), replace = TRUE)
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
	ggplot(aes(x = vb, y = LL, group = n)) + 
	geom_line(stat = "identity", alpha = .35, color = "#333333", size = .25) +
	theme_classic() +
	xlab(bquote(atop(" ", nu[b]))) +
	ylab("\nSymmetric Log Likelihood\n") +
	scale_x_continuous(breaks = 1:7,
			   labels = 1:7) +
	scale_y_continuous(labels = scientific_10)

pdf(file = "vb_LL.pdf", height = 3, width = 3)
print(plot_)
dev.off()

LL = vector(mode = "numeric", length = 7)
for (i in 1:7) {
	LL[i] = SymmLL(m = m, c = c, nb = i)$LL
}

P = 1 - pchisq(2*(LL[2]-LL[1]),1)
P = 1 - pchisq(2*(LL[3]-LL[2]),1)
P = 1 - pchisq(2*(LL[4]-LL[3]),1)
P = 1 - pchisq(2*(LL[5]-LL[4]),1)

LL = SymmLL(m = m, c = c, nb = 5)

data_ = do.call(rbind, LL$p_bjr) %>%
	dplyr::as_tibble() %>%
	dplyr::mutate(UUID = rep(paste(data$Case_ID, ":", data$Gene_Symbol, ":", data$HGVSp_Short), .MMEnv$n_run),
		      VAF = rep(data$N_Alt/data$N_Total, .MMEnv$n_run)) %>%
	dplyr::arrange(VAF) %>%
	dplyr::mutate(UUID = factor(UUID, levels = unique(UUID), ordered = TRUE))

plot_ = data_ %>%
	ggplot(aes(x = V1, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = hex_cols[1], color = hex_cols[1], alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==1)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

pdf(file = "p(vb=1).pdf", height = 10, width = 5)
print(plot_)
dev.off()

plot_ = data_ %>%
	ggplot(aes(x = V2, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = hex_cols[2], color = hex_cols[2], alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==2)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

pdf(file = "p(vb=2).pdf", height = 10, width = 5)
print(plot_)
dev.off()

plot_ = data_ %>%
	ggplot(aes(x = V3, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = hex_cols[3], color = hex_cols[3], alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==3)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

pdf(file = "p(vb=3).pdf", height = 10, width = 5)
print(plot_)
dev.off()

plot_ = data_ %>%
	ggplot(aes(x = V4, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = hex_cols[4], color = hex_cols[4], alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==4)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

pdf(file = "p(vb=4).pdf", height = 10, width = 5)
print(plot_)
dev.off()

plot_ = data_ %>%
	ggplot(aes(x = V5, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = hex_cols[5], color = hex_cols[5], alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==5)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

pdf(file = "p(vb=5).pdf", height = 10, width = 5)
print(plot_)
dev.off()

data_ = dplyr::tibble(vb = apply(LL$p_bj, 1, which.max)) %>%
	dplyr::mutate(UUID = paste0(data$Gene_Symbol, " ", gsub("p.", "", data$HGVSp_Short)),
		      VAF = data$N_Alt/data$N_Total) %>%
	dplyr::arrange(desc(VAF)) %>%
	dplyr::mutate(UUID = factor(UUID, levels = unique(UUID), ordered = TRUE))

plot_ = data_ %>%
	ggplot(aes(x = UUID, y = 100*VAF, fill = factor(vb))) +
	geom_hline(yintercept = .MMEnv$vb[1:5]*100, colour = "black", linetype = 1, size = .25) +
	geom_point(stat = "identity", shape = 21, size = 2) +
	scale_fill_manual(values = hex_cols) +
	xlab("\n\n") +
	ylab("\nVAF (%)\n\n") +
	scale_y_continuous(breaks = seq(from = 0, to = 25, by = 5),
			   labels = seq(from = 0, to = 25, by = 5)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
	guides(fill = FALSE)
	
pdf(file = "VAF_by_Variant.pdf", width = 6, height = 5)
print(plot_)
dev.off()
