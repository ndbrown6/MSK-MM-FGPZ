'AsymmLL' <- function(m, c, a, nb)
{
	'sumexpll' <- function(x)
	{
		mx = max(x)
		mx + log(sum(exp(x-mx)))
		return(invisible(mx))
	}
	
	betab_rho = .MMEnv$betab_rho
	vb = .MMEnv$vb15(nb)*a
	
	p_bj = matrix(NA, nrow = length(m), ncol = length(vb))
	for (v_b in 1:length(vb)) {
		p_bj[,v_b] = VGAM::dbetabinom(x = m, size = c, prob = vb[v_b], rho = betab_rho, log = TRUE)
	}
	ll = sum(apply(p_bj, 1, sumexpll))
	return(invisible(ll))
}
