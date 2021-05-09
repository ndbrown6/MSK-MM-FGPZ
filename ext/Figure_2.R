#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
library('mosaicm')
library('ggplot2')
library('ggridges')

# test symmetric model with actual data vb = n and nb = 4 branches
data("vb=n")
m = data %>% .[["N_Alt"]]
c = data %>% .[["N_Total"]]
LL = SymmLL(m = m, c = c, nb = 4)

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

plot_ = data_ %>%
	ggplot(aes(x = V2, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = "#fdae61", color = "#fdae61", alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==2)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

plot_ = data_ %>%
	ggplot(aes(x = V3, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = "#008837", color = "#008837", alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==3)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

plot_ = data_ %>%
	ggplot(aes(x = V4, y = UUID, group = UUID)) + 
	geom_density_ridges(stat = "density_ridges", fill = "#2c7bb6", color = "#2c7bb6", alpha = .75) +
	theme_classic() +
	xlab(bquote(atop(" ", Pr(nu[b] ==4)))) +
	ylab("") +
	scale_x_continuous(limits = c(-0.1,1.1),
			   breaks = c(0, .2, .4, .6, .8, 1))

# test symmetric model with actual data vb = n and nb = 15 branches
data("vb=n")
m = data %>% .[["N_Alt"]]
c = data %>% .[["N_Total"]]
LL = SymmLL(m = m, c = c, nb = 15)

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
a1 = seq(from = .1, to = 2, length = 100)
a2 = seq(from = .1, to = 4, length = 100)
a3 = seq(from = .1, to = 8, length = 100)
a4 = seq(from = .1, to = 16, length = 100)

LLvb1 = LLvb2 = LLvb3 = LLvb4 = matrix(NA, nrow = 100, ncol = 100)

pb = txtProgressBar(min = 1, max = 100, style = 3)
for (i in 1:100) {
	setTxtProgressBar(pb, i)
	for (j in 1:100) {
		LLvb1[i, j] = AsymmLL(m = m, c = c, a = c(a1[i], rep(1,15)), nb = 15)$LL -
			      AsymmLL(m = m, c = c, a = c(1, a2[j], rep(1,13)), nb = 15)$LL
		
		LLvb2[i, j] = AsymmLL(m = m, c = c, a = c(1, a2[i], rep(1,13)), nb = 15)$LL -
			      AsymmLL(m = m, c = c, a = c(rep(1,3), a3[j], rep(1,11)), nb = 15)$LL
		
		LLvb3[i, j] = AsymmLL(m = m, c = c, a = c(1, a2[i], rep(1,13)), nb = 15)$LL -
			      AsymmLL(m = m, c = c, a = c(rep(1,7), a4[j], rep(1,7)), nb = 15)$LL
		
		LLvb4[i, j] = AsymmLL(m = m, c = c, a = c(rep(1,3), a3[i], rep(1,11)), nb = 15)$LL -
			      AsymmLL(m = m, c = c, a = c(rep(1,7), a4[j], rep(1,7)), nb = 15)$LL
	}
}
close(pb)
