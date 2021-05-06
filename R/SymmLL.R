'SymmLL' <- function(m_j, c_j, b, method = "binomial", log = TRUE)
{
	if (log) {
		if (method == "binomial") {
			ll = dbinom(x = m_j, size = c_j, prob = .MMEnv$v_b[b], log = TRUE) +
		     	     .MMEnv$n_b[b]
		} else if (method == "betabinomial") {
			ll = dbeta(x=m_j/c_j, round(c_j*.MMEnv$v_b[b])+1, round(c_j*(1-.MMEnv$v_b[b]))+1, log = TRUE) +
			     .MMEnv$n_b[b]
		}
	} else {
		if (method == "binomial") {
			ll = dbinom(x = m_j, size = c_j, prob = .MMEnv$v_b[b]) *
		     	     .MMEnv$n_b[b]
		} else if (method == "betabinomial") {
			ll = dbeta(x=m_j/c_j, round(c_j*.MMEnv$v_b[b])+1, round(c_j*(1-.MMEnv$v_b[b]))+1) *
			     .MMEnv$n_b[b]
		}
	}
	return(invisible(ll))
}
