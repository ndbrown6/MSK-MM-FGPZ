'LL' <- function(m_j, c_j, n, b) {
	ll = 0
	for (j in 1:n) {
		ll = ll + dbinom(x = m_j[j], size = c_j[j], prob = .MMEnv$v_b[b], log = TRUE) + dbinom(x = c_j[j] - m_j[j], size = c_j[j], prob = 1 - .MMEnv$v_b[b], log = TRUE)
	}
	return(invisible(ll))
}
