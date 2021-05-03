'AsymmLL' <- function(m_j, c_j, b, a, method = "binomial", log = TRUE)
{
	if (log) {
		if (method == "binomial") {
			ll = dbinom(x = m_j, size = c_j, prob = .MMEnv$v_b[b]*a, log = TRUE) +
			     dbinom(x = c_j - m_j, size = c_j, prob = (1 - .MMEnv$v_b[b]*a), log = TRUE) +
		     	     .MMEnv$n_b[b]
		} else if (method == "betabinomial") {
			ll = dbeta(x=m_j/c_j, round(c_j*.MMEnv$v_b[b]*a)+1, round(c_j*(1-(.MMEnv$v_b[b]*a)))+1, log = TRUE) +
			     dbeta(x=(c_j-m_j)/c_j, round(c_j*(1-(.MMEnv$v_b[b]*a)))+1, round(c_j*.MMEnv$v_b[b]*a)+1, log = TRUE) +
			     .MMEnv$n_b[b]
		}
	} else {
		if (method == "binomial") {
			ll = dbinom(x = m_j, size = c_j, prob = .MMEnv$v_b[b]*a) *
			     dbinom(x = c_j - m_j, size = c_j, prob = 1 - (.MMEnv$v_b[b]*a)) *
		     	     .MMEnv$n_b[b]
		} else if (method == "betabinomial") {
			ll = dbeta(x=m_j/c_j, round(c_j*.MMEnv$v_b[b]*a)+1, round(c_j*(1-(.MMEnv$v_b[b]*a)))+1) *
			     dbeta(x=(c_j-m_j)/c_j, round(c_j*(1-(.MMEnv$v_b[b]*a)))+1, round(c_j*.MMEnv$v_b[b]*a)+1) *
			     .MMEnv$n_b[b]
		}
	}
	return(invisible(ll))
}
