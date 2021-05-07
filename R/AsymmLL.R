'AsymmLL' <- function(m_j, c_j, v_b, a)
{
	ll = VGAM::dbetabinom(x = m_j, size = c_j, prob = v_b*a, rho = 1e-03)
	return(invisible(ll))
}
