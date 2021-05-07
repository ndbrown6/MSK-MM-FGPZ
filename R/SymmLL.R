'SymmLL' <- function(m_j, c_j, v_b)
{
	ll = VGAM::dbetabinom(x = m_j, size = c_j, prob = v_b, rho = 1e-04)
	return(invisible(ll))
}
