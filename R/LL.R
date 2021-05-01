'LL' <- function(m_j, c_j, b) {
	ll = dbinom(x = m_j, size = c_j, prob = .MMEnv$v_b[b]) * dbinom(x = c_j - m_j, size = c_j, prob = .MMEnv$v_b[b])
	return(invisible(ll))
}

