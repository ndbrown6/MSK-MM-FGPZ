#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
library('mskmmfgpz')
library('reshape2')

'scientific_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^", scales::scientific_format()(x)), , fixed = TRUE))
}


data("vb=1")

ll = list()
ll[[1]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[1]][i] = Total_LL(m, c, n, i, method = "binomial", log = TRUE)
}

ll[[2]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[2]][i] = Total_LL(m, c, n, i, method = "betabinomial", log = TRUE)
}

ll[[3]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[3]][i] = Total_LL(m, c, n, i, method = "binomial", log = FALSE)
}

ll[[4]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[4]][i] = Total_LL(m, c, n, i, method = "betabinomial", log = FALSE)
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
	facet_wrap(~method+log, scales = "free_y")

pdf(file = "vb=1.pdf", width = 7, height = 5)
print(plot_)
dev.off()


data("vb=2")

ll = list()
ll[[1]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[1]][i] = Total_LL(m, c, n, i, method = "binomial", log = TRUE)
}

ll[[2]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[2]][i] = Total_LL(m, c, n, i, method = "betabinomial", log = TRUE)
}

ll[[3]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[3]][i] = Total_LL(m, c, n, i, method = "binomial", log = FALSE)
}

ll[[4]] = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll[[4]][i] = Total_LL(m, c, n, i, method = "betabinomial", log = FALSE)
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
	facet_wrap(~method+log, scales = "free_y")

pdf(file = "vb=2.pdf", width = 7, height = 5)
print(plot_)
dev.off()

'A_LL' <- function(m_j, c_j, b, a_1)
{
	ll = dbeta(x=m_j/c_j, round(c_j*.MMEnv$v_b[b])+1, round(c_j*(1-.MMEnv$v_b[b]))+1, log = TRUE) +
	     dbeta(x=(c_j-m_j)/c_j, round(c_j*(1-.MMEnv$v_b[b]))+1, round(c_j*.MMEnv$v_b[b])+1, log = TRUE) +
	     .MMEnv$n_b[b]
	return(invisible(ll))
}
