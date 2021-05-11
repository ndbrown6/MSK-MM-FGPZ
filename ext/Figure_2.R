#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
library('mosaicm')
library('ggplot2')
library('ggridges')

data("vb=n")
m = data %>% .[["N_Alt"]]
c = data %>% .[["N_Total"]]
hex_cols = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")

# test symmetric model with actual data vb = n and nb = 3|4|5|6 cell generations
LL3 = SymmLL(m = m, c = c, nb = 3)
LL4 = SymmLL(m = m, c = c, nb = 4)
LL5 = SymmLL(m = m, c = c, nb = 5)
LL6 = SymmLL(m = m, c = c, nb = 6)
P = 1 - pchisq(2*(LL5$LL-LL4$LL),1)

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


# test asymmetric model with actual data vb = n and nb = 5 cell generations
n = 25
a = seq(from = 0, to = 1, length = n)
LL = matrix(NA, nrow = n, ncol = n)

pb = txtProgressBar(min = 1, max = n, style = 3)
for (i in 1:n) {
	setTxtProgressBar(pb, i)
	for (j in 1:n) {
		ai = rep(1, 15)
		ai[4] = a[i]
		ai[8] = a[j]
		LL[i, j] = AsymmLL(m = m, c = c, a = ai)$LL
	}
}
close(pb)
