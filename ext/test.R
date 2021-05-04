#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
library('mosaicm')
library('ggforce')

'scientific_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^", scales::scientific_format()(x)), , fixed = TRUE))
}


# test functions with simulated data at v_b = 1
data("vb=1")

ll = list()
ll[[1]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[1]][i] = TotalLL(m, c, i, method = "binomial", log = TRUE)
}

ll[[2]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[2]][i] = TotalLL(m, c, i, method = "betabinomial", log = TRUE)
}

ll[[3]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[3]][i] = TotalLL(m, c, i, method = "binomial", log = FALSE)
}

ll[[4]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[4]][i] = TotalLL(m, c, i, method = "betabinomial", log = FALSE)
}

ll = do.call(cbind, ll) %>%
     reshape2::melt() %>%
     dplyr::as_tibble() %>%
     dplyr::select(LL = value) %>%
     dplyr::mutate(b = rep(1:4, times = 4),
		   method = rep(c("Binomial", "Beta-Binomial"), each = 4, times = 2),
		   log = rep(c("Log", "Non-Log"), each = 8))

plot_ = ll %>%
	dplyr::mutate(LL = ifelse(LL>1E50, 1E50, LL)) %>%
	dplyr::mutate(LL = ifelse(LL<(-1E50), -1E50, LL)) %>%
	ggplot(aes(x = b, y = LL)) +
	geom_point(stat = "identity", shape = 21, fill = "white", size = 2) +
	xlab("\n\nb\n") +
	ylab("\nLL\n\n") +
	scale_y_continuous(labels = scientific_10) +
	facet_wrap(method~log, scales = "free_y")

pdf(file = "vb=1.pdf", width = 7, height = 5)
print(plot_)
dev.off()

# test functions with simulated data at v_b = 2
data("vb=2")

ll = list()
ll[[1]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[1]][i] = TotalLL(m, c, i, method = "binomial", log = TRUE)
}

ll[[2]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[2]][i] = TotalLL(m, c, i, method = "betabinomial", log = TRUE)
}

ll[[3]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[3]][i] = TotalLL(m, c, i, method = "binomial", log = FALSE)
}

ll[[4]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[4]][i] = TotalLL(m, c, i, method = "betabinomial", log = FALSE)
}

ll = do.call(cbind, ll) %>%
     reshape2::melt() %>%
     dplyr::as_tibble() %>%
     dplyr::select(LL = value) %>%
     dplyr::mutate(b = rep(1:4, times = 4),
		   method = rep(c("Binomial", "Beta-Binomial"), each = 4, times = 2),
		   log = rep(c("Log", "Non-Log"), each = 8))

plot_ = ll %>%
	dplyr::mutate(LL = ifelse(LL>1E50, 1E50, LL)) %>%
	dplyr::mutate(LL = ifelse(LL<(-1E50), -1E50, LL)) %>%
	ggplot(aes(x = b, y = LL)) +
	geom_point(stat = "identity", shape = 21, fill = "white", size = 2) +
	xlab("\n\nb\n") +
	ylab("\nLL\n\n") +
	scale_y_continuous(labels = scientific_10) +
	facet_wrap(method~log, scales = "free_y")

pdf(file = "vb=2.pdf", width = 7, height = 5)
print(plot_)
dev.off()

# test asymmetric cell divisions with simulated data at v_b = 1
data("vb=1")
eps = c(1, 5, 10, 20, 30, 50)
ll = list()
for (j in 1:length(eps)) {
	mz = m + eps[j]
	mz[(mz/c) > 1] = c[(mz/c) > 1]
	a = seq(from = 0.1, to = 2, length = 100)
	ll[[j]] = vector(mode = "numeric", length = 100)
	for (i in 1:100) {
		ll[[j]][i] = TotalLL(mz, c, 1, a = a[i], method = "betabinomial", log = TRUE)
	}
}
ll = dplyr::tibble(LL = unlist(ll),
		   alpha = rep(seq(from = 0.1, to = 2, length = 100), times = length(eps)),
		   eps = rep(eps, each = 100))

plot_ = ll %>%
	ggplot(aes(x = alpha, y = LL, color = factor(eps), group = eps)) +
	geom_line(stat = "identity", size = 1, alpha = .55) +
	xlab(bquote(atop(" ", alpha))) +
	ylab("\nLL\n\n") +
	scale_y_continuous(labels = scientific_10) +
	facet_zoom(xlim = c(.75, 2), ylim = c(-10000, 0), horizontal = TRUE, show.area = TRUE) +
	guides(color = guide_legend(title = expression(epsilon)))

pdf(file = "vb=1_eps.pdf", width = 10, height = 5)
print(plot_)
dev.off()

# test asymmetric cell divisions with simulated data at v_b = 2
data("vb=2")
eps = c(1, 5, 10, 20, 30, 50)
ll = list()
for (j in 1:length(eps)) {
	mz = m + eps[j]
	mz[(mz/c) > 1] = c[(mz/c) > 1]
	a = seq(from = 0.1, to = 4, length = 100)
	ll[[j]] = vector(mode = "numeric", length = 100)
	for (i in 1:100) {
		ll[[j]][i] = TotalLL(mz, c, 2, a = a[i], method = "betabinomial", log = TRUE)
	}
}
ll = dplyr::tibble(LL = unlist(ll),
		   alpha = rep(seq(from = 0.1, to = 4, length = 100), times = length(eps)),
		   eps = rep(eps, each = 100))

plot_ = ll %>%
	ggplot(aes(x = alpha, y = LL, color = factor(eps), group = eps)) +
	geom_line(stat = "identity", size = 1, alpha = .55) +
	xlab(bquote(atop(" ", alpha))) +
	ylab("\nLL\n\n") +
	scale_y_continuous(labels = scientific_10) +
	facet_zoom(xlim = c(.75, 3), ylim = c(-10000, 0), horizontal = TRUE, show.area = TRUE) +
	guides(color = guide_legend(title = expression(epsilon)))

pdf(file = "vb=2_eps.pdf", width = 10, height = 5)
print(plot_)
dev.off()
