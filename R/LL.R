'LL' <- function(m_j, c_j, b, alpha = .5) {
	ll = alpha*dbinom(x = m_j, size = c_j, prob = .MMEnv$v_b[b]) * (1-alpha)*dbinom(x = c_j - m_j, size = c_j, prob = 1 - .MMEnv$v_b[b])
	return(invisible(ll))
}

