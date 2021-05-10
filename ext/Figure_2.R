#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
library('mosaicm')
library('ggplot2')
library('ggridges')

# test symmetric model with actual data vb = n and nb = 15 branches
data("vb=n")
m = data %>% .[["N_Alt"]]
c = data %>% .[["N_Total"]]
LL = SymmLL(m = m, c = c)

data_ = do.call(rbind, LL$p_bjr) %>%
	dplyr::as_tibble() %>%
	dplyr::mutate(UUID = rep(paste(data$Case_ID, ":", data$Gene_Symbol, ":", data$HGVSp_Short), .MMEnv$n_run),
		      VAF = rep(data$N_Alt/data$N_Total, .MMEnv$n_run)) %>%
	dplyr::arrange(VAF) %>%
	dplyr::mutate(UUID = factor(UUID, levels = unique(UUID), ordered = TRUE))

plot_ = data_ %>%
	ggplot(aes(x = V1, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = "#d7191c", color = "#d7191c", alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==1)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

pdf(file = "p(vb=1).pdf", height = 10, width = 5)
print(plot_)
dev.off()

plot_ = data_ %>%
	ggplot(aes(x = V2+V3, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = "#fdae61", color = "#fdae61", alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==2)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

pdf(file = "p(vb=2).pdf", height = 10, width = 5)
print(plot_)
dev.off()

plot_ = data_ %>%
	ggplot(aes(x = V4+V5+V6+V7, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = "#008837", color = "#008837", alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==3)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

pdf(file = "p(vb=3).pdf", height = 10, width = 5)
print(plot_)
dev.off()

plot_ = data_ %>%
	ggplot(aes(x = V8+V9+V10+V11+V12+V13+V14+V15, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = "#2c7bb6", color = "#2c7bb6", alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==4)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

pdf(file = "p(vb=4).pdf", height = 10, width = 5)
print(plot_)
dev.off()

# test asymmetric model with actual data vb = n and nb = 15 branches
n = 25
a1 = seq(from = 0, to = 2, length = n)
a2 = seq(from = 0, to = 4, length = n)
LL = matrix(NA, nrow = n, ncol = n)

pb = txtProgressBar(min = 1, max = n, style = 3)
for (i in 1:n) {
	setTxtProgressBar(pb, i)
	for (j in 1:n) {
		LL[i, j] = AsymmLL(m = m, c = c, a = c(a1[i], 1, rep(1,13)))$LL - AsymmLL(m = m, c = c, a = c(1, a2[j], rep(1,13)))$LL
	}
}
close(pb)
