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

'fsq' <- function(alpha, beta, s_q, q_t)
{
	return(invisible(((alpha*s_q) + beta)/((alpha*q_t)+(1-alpha)*2)))
}

alpha = seq(from = 0, to = .5, by = .01)
s_q = c(1, 1, 2, 1, 2, 3)
q_t = c(1, 2, 2, 3, 3, 3)
col = c("steelblue", "goldenrod3", "salmon")
lty = c(1, 2, 3)

beta = 0.5
data_ = list()
for (i in 1:length(s_q)) {
	data_[[i]] = dplyr::tibble(alpha = alpha,
				   s_q = s_q[i],
				   q_t = q_t[i],
				   f_sq = fsq(alpha = alpha, beta = beta, s_q = s_q[i], q_t = q_t[i]),
				   beta = paste0("Beta = ", beta))
}
data_ = do.call(rbind, data_)

plot_ = data_ %>%
	ggplot(aes(x = alpha, y = f_sq, color = factor(q_t), linetype = factor(s_q))) +
	geom_line(stat = "identity", size = 1) +
	geom_vline(xintercept = .2, linetype = 3) +
	geom_hline(yintercept = fsq(alpha = .2, beta = beta, s_q = 1, q_t = 3), linetype = 3) +
	xlab(expression(alpha)) +
	ylab(expression(f[sq])) +
	scale_x_continuous(limits = c(0,.5)) +
	scale_y_continuous(limits = c(0,1)) +
	guides(color = guide_legend(title = expression("Total copies (q"[t]~")"), order = 1),
	       linetype = guide_legend(title = expression("Mutated copies (s"[q]~")"), order = 2)) +
	facet_wrap(~beta)

pdf(file = "beta=.5.pdf", width = 5, height = 3)
print(plot_)
dev.off()

beta = 0.25
data_ = list()
for (i in 1:length(s_q)) {
	data_[[i]] = dplyr::tibble(alpha = alpha,
				   s_q = s_q[i],
				   q_t = q_t[i],
				   f_sq = fsq(alpha = alpha, beta = beta, s_q = s_q[i], q_t = q_t[i]),
				   beta = paste0("Beta = ", beta))
}
data_ = do.call(rbind, data_)

plot_ = data_ %>%
	ggplot(aes(x = alpha, y = f_sq, color = factor(q_t), linetype = factor(s_q))) +
	geom_line(stat = "identity", size = 1) +
	geom_vline(xintercept = .2, linetype = 3) +
	geom_hline(yintercept = fsq(alpha = .2, beta = beta, s_q = 1, q_t = 3), linetype = 3) +
	xlab(expression(alpha)) +
	ylab(expression(f[sq])) +
	scale_x_continuous(limits = c(0,.5)) +
	scale_y_continuous(limits = c(0,1)) +
	guides(color = guide_legend(title = expression("Total copies (q"[t]~")"), order = 1),
	       linetype = guide_legend(title = expression("Mutated copies (s"[q]~")"), order = 2)) +
	facet_wrap(~beta)


pdf(file = "beta=.25.pdf", width = 5, height = 3)
print(plot_)
dev.off()

beta = 0.125
data_ = list()
for (i in 1:length(s_q)) {
	data_[[i]] = dplyr::tibble(alpha = alpha,
				   s_q = s_q[i],
				   q_t = q_t[i],
				   f_sq = fsq(alpha = alpha, beta = beta, s_q = s_q[i], q_t = q_t[i]),
				   beta = paste0("Beta = ", beta))
}
data_ = do.call(rbind, data_)

plot_ = data_ %>%
	ggplot(aes(x = alpha, y = f_sq, color = factor(q_t), linetype = factor(s_q))) +
	geom_line(stat = "identity", size = 1) +
	geom_vline(xintercept = .2, linetype = 3) +
	geom_hline(yintercept = fsq(alpha = .2, beta = beta, s_q = 1, q_t = 3), linetype = 3) +
	xlab(expression(alpha)) +
	ylab(expression(f[sq])) +
	scale_x_continuous(limits = c(0,.5)) +
	scale_y_continuous(limits = c(0,1)) +
	guides(color = guide_legend(title = expression("Total copies (q"[t]~")"), order = 1),
	       linetype = guide_legend(title = expression("Mutated copies (s"[q]~")"), order = 2)) +
	facet_wrap(~beta)

pdf(file = "beta=.125.pdf", width = 5, height = 3)
print(plot_)
dev.off()


alpha = seq(from = 0, to = .5, by = .01)
beta = seq(from = 0.01, to = .6, by = .1)
s_q = 1
q_t = 3
data_ = list()
for (i in 1:length(beta)) {
	data_[[i]] = dplyr::tibble(alpha = alpha,
				   beta = beta[i],
				   f_sq = fsq(alpha = alpha, beta = beta[i], s_q = s_q, q_t = q_t))
		
}
data_ = do.call(rbind, data_)

plot_ = data_ %>%
	ggplot(aes(x = alpha, y = f_sq, group = factor(beta))) +
	geom_line(stat = "identity", size = 1, color = "#772953") +
	geom_vline(xintercept = .2, linetype = 3) +
	geom_hline(yintercept = fsq(alpha = .2, beta = .01, s_q = s_q, q_t = q_t), linetype = 3) +
	xlab(expression(alpha)) +
	ylab(expression(f[sq])) +
	scale_x_continuous(limits = c(0,.5)) +
	scale_y_continuous(limits = c(0,.5)) +
	annotate(geom = "rect", xmin = .355, xmax = .395, ymin = .15, ymax = .4, fill = "gray92", alpha = .99) +
	geom_text(x = .3775, y = .16, label = "0.01", size = 2) +
	geom_text(x = .3775, y = .2, label = "0.10", size = 2) +
	geom_text(x = .3775, y = .25, label = "0.20", size = 2) +
	geom_text(x = .3775, y = .29, label = "0.30", size = 2) +
	geom_text(x = .3775, y = .335, label = "0.40", size = 2) +
	geom_text(x = .3775, y = .375, label = "0.50", size = 2)

pdf(file = "sq=1_qt=3.pdf", width = 3.25, height = 3)
print(plot_)
dev.off()
