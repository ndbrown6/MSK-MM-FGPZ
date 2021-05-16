'SymmLL' <- function(m, c, nb)
{
	n_run = .MMEnv$n_run
	betab_rho = .MMEnv$betab_rho
	max_iter = .MMEnv$max_iter
	min_dist = .MMEnv$min_dist
	vb = .MMEnv$vb[1:nb]
	
	LL_b = -Inf
	p_bjr = r_br = list()
	for (r in 1:n_run) {
		rb_old = rep(-1, times = nb)
		rb_est = rep(length(m)/nb, times = nb)
		iter = 0
		while (iter<max_iter & sum(abs(rb_est-rb_old))>min_dist) {
			iter = iter + 1
			rb_old = rb_est
			p_bj = matrix(NA, nrow = length(m), ncol = nb)
			for (v_b in 1:nb) {
				p_bj[,v_b] = VGAM::dbetabinom(x = m, size = c, prob = vb[v_b], rho = betab_rho) * rb_est[v_b] * 2^v_b
			}
			p_bj = pmax(p_bj, 1e-20)
			p_bj = p_bj/apply(p_bj, 1, sum)
			rb_est = apply(p_bj, 2, sum)
		}
		rb = apply(p_bj, 2, sum)
		loglik = matrix(NA, nrow = length(m), ncol = nb)
		for (v_b in 1:nb) {
			loglik[,v_b] = log(rb[v_b]) + log(2^v_b) + VGAM::dbetabinom(x = m, size = c, prob = vb[v_b], rho = betab_rho, log = TRUE)
		}
		LL = sum(apply(loglik, 1, function(x) { log(sum(exp(x))/sum(rb_est)) }))
		if (LL>LL_b) {
			LL_b = LL
			params = list(LL = LL,
				      r_b = rb,
				      p_bj = p_bj)
		}
		p_bjr[[r]] = p_bj
		r_br[[r]] = rb_est
	}
	params = list(LL = params$LL,
		      r_b = params$r_b,
		      p_bj = params$p_bj,
		      r_br = r_br,
		      p_bjr = p_bjr)
	return(invisible(params))

}
