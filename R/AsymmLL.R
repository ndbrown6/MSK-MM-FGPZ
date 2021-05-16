'AsymmLL' <- function(m, c, a, nb)
{
	betab_rho = .MMEnv$betab_rho
	vb = .MMEnv$vb15(nb)*a
	
	p_bj = matrix(NA, nrow = length(m), ncol = length(vb))
	for (v_b in 1:length(vb)) {
		p_bj[,v_b] = VGAM::dbetabinom(x = m, size = c, prob = vb[v_b], rho = betab_rho, log = TRUE)
	}
	
	'logsumexp' <- function(llvec)
	{
		maxllvec = max(llvec)
		maxllvec + log(sum(exp(llvec-maxllvec)))
		return(invisible(maxllvec))
	}
	
	return(invisible(sum(apply(p_bj, 1, logsumexp))))

}
