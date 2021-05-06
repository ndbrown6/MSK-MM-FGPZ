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
data("vb=2")1, method = "betabinomial", log = TRUE - t v_b = 2
data("vb=2")1, method = "betabinomial", log = TRUE
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


##################################################################
'SymmLL' <- function(m_j, c_j, v_b)
{
	ll = dbinom(x = m_j, size = c_j, prob = v_b, log = TRUE)
	return(invisible(ll))
}

'TotalLL' <- function(m, c, v)
{
	ll = sum(AsymmLL(m, c, v))
	return(invisible(ll))
}



# test symmetric functions with simulated data at v_b = 1
data("vb=1")

ll = vector(mode = "numeric", length = 15)
vb = c(rep(.25, 1), rep(.125, 2), rep(.0625, 4), rep(.03125, 8))
for (i in 1:15) {
	ll[i] = TotalLL(m, c, vb[i])
}

# test symmetric functions with simulated data at v_b = 2
data("vb=2")

ll = vector(mode = "numeric", length = 15)
vb = c(rep(.25, 1), rep(.125, 2), rep(.0625, 4), rep(.03125, 8))
for (i in 1:15) {
	ll[i] = TotalLL(m, c, vb[i])
}

##################################################################
'AsymmLL' <- function(m_j, c_j, v_b, a)
{
	ll = dbinom(x = m_j, size = c_j, prob = v_b, log = TRUE) + a
	return(invisible(ll))
}

'TotalLL' <- function(m, c, v, a)
{
	ll = sum(AsymmLL(m, c, v, a))
	return(invisible(ll))
}

# test asymmetric functions with simulated data
data("vb=1")
l.m = m[1:25]
l.c = c[1:25]

data("vb=2")
m = c(l.m, m[1:25])
c = c(l.c, c[1:25])

n = 50

load("~/GitHub/MSK-MM-FGPZ/data/vb=0.RData")
m = data_$m
c = data_$c
n = 38

vb = c(rep(.25, 1), rep(.125, 2), rep(.0625, 4), rep(.03125, 8))
LL = matrix(NA, nrow = 100, ncol = 100)
a1 = seq(from = 1, to = 2, length = 100)
a2 = seq(from = 1, to = 4, length = 100)
for (ii in 1:100) {
	a_1 = c(a1[ii], rep(1, 14))
	for (jj in 1:100) {
		a_2 = c(1, a2[jj], rep(1, 13))
		ll0 = ll1 = vector(mode = "numeric", length = n)
		for (i in 1:n) {
			for (j in 1:15) {
				ll0[i] = ll0[i] + TotalLL(m[i], c[i], vb[j], a_1[j])
				ll1[i] = ll1[i] + TotalLL(m[i], c[i], vb[j], a_2[j])
			}
		}
		LL[ii, jj] = sum(ll0) - sum(ll1)
	}
}
