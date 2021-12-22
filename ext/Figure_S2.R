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
library("pwr")

'scientific_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^", scales::scientific_format()(x)), , fixed = TRUE))
}


# vaf versus coverage
n_len = 1000
p1 = seq(from = 0.01125, to = 0.25, length = n_len)
p2 = 10E-4
power = c(0.8, 0.9, 0.95, 0.99, 0.999)
n  = list()

for (i in 1:length(power)) {
	
	n[[i]] = vector(mode = "numeric", length = n_len)
	
	for (j in 1:n_len) {
		
		n[[i]][j] = pwr.p.test(h = ES.h(p1 = p1[j], p2 = p2), sig.level = 0.05, power = power[i], alternative = "greater")$n
	}
	
}
data_ = do.call(cbind, n) %>%
	dplyr::as_tibble() %>%
	dplyr::rename(`80.0` = V1,
		      `90.0` = V2,
		      `95.0` = V3,
		      `99.0` = V4,
		      `99.9` = V5) %>%
	dplyr::mutate(VAF = p1*100) %>%
	reshape2::melt(id.vars = "VAF", value.name = "Coverage") %>%
	dplyr::as_tibble()

plot_ = data_ %>%
	ggplot(aes(x = VAF, y = Coverage, color = variable)) +
	geom_line(stat = "identity", size = 1) +
	xlab("\n\nVAF (%)\n") +
	ylab("\nCoverage\n\n") +
	scale_x_continuous() +
	scale_y_log10(labels = scientific_10) +
	scale_color_brewer(type = "seq", palette = "GnBu") +
	theme_classic() +
	guides(color = guide_legend(title = "Detection\nProbability (%)")) +
	annotation_logticks(side = "l")

pdf(file = "nbyvaf.pdf", width = 6, height = 5)
print(plot_)
dev.off()

# vaf versus power
n_len = 1000
p1 = seq(from = 0.001, to = 0.25, length = n_len)
p2 = 10E-4
n = c(20, 30, 40, 50, 100, 500, 1000)
power  = list()

for (i in 1:length(n)) {
	
	power[[i]] = vector(mode = "numeric", length = n_len)
	
	for (j in 1:n_len) {
		
		power[[i]][j] = pwr.p.test(h = ES.h(p1 = p1[j], p2 = p2), sig.level = 0.05, n = n[i], alternative = "greater")$power
	}
	
}
data_ = do.call(cbind, power) %>%
	dplyr::as_tibble() %>%
	dplyr::rename(`20` = V1,
		      `30` = V2,
		      `40` = V3,
		      `50` = V4,
		      `100` = V5,
		      `500` = V6,
		      `1000` = V7) %>%
	dplyr::mutate(VAF = p1*100) %>%
	reshape2::melt(id.vars = "VAF", value.name = "Probability") %>%
	dplyr::as_tibble()

plot_ = data_ %>%
	ggplot(aes(x = VAF, y = Probability, color = variable)) +
	geom_line(stat = "identity", size = 1) +
	xlab("\n\nVAF (%)\n") +
	ylab("\nProbability\n\n") +
	scale_x_log10() +
	scale_y_continuous() +
	scale_color_brewer(type = "seq", palette = "GnBu") +
	theme_classic() +
	guides(color = guide_legend(title = "Sequence      \nDepth")) +
	annotation_logticks(side = "b")

pdf(file = "powerbyvaf.pdf", width = 6, height = 5)
print(plot_)
dev.off()

# coverage versus power
n_len = 1000
p1 = 0.015
p2 = 10E-4
power = .8

n = pwr.p.test(h = ES.h(p1 = p1, p2 = p2), sig.level = 0.05, power = power, alternative = "greater")$n
power = pwr.p.test(h = ES.h(p1 = p1, p2 = p2), sig.level = 0.05, n = 250, alternative = "greater")$power