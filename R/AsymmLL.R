'AsymmLL' <- function(m, c, a)
{
	n_run = .MMEnv$n_run
	betab_rho = .MMEnv$betab_rho
	max_iter = .MMEnv$max_iter
	min_dist = .MMEnv$min_dist
	nb = .MMEnv$nb
	vb = .MMEnv$vb*a
	
	LL_b = -Inf
	p_bjr = list()
	for (r in 1:n_run) {
		rb_old = rep(-1, times = nb)
		rb_est = rep(1/nb, times = nb)
		iter = 0
		while (iter<max_iter & sum(abs(rb_est-rb_old))>min_dist) {
			iter = iter + 1
			rb_old = rb_est
			p_bj = matrix(NA, nrow = length(m), ncol = nb)
			for (v_b in 1:nb) {
				p_bj[,v_b] = VGAM::dbetabinom(x = m, size = c, prob = vb[v_b], rho = betab_rho) * rb_est[v_b]
			}
			p_bj = pmax(p_bj, 1e-20)
			p_bj = p_bj/apply(p_bj, 1, sum)
			for (v_b in 1:nb) {
				rb_est[v_b] = sum(m*p_bj[,v_b])/sum(c*p_bj[,v_b])
			}
		}
		rb = apply(p_bj, 2, sum)/length(m)
		loglik = matrix(NA, nrow = length(m), ncol = nb)
		for (v_b in 1:nb) {
			loglik[,v_b] = log(rb[v_b]) + VGAM::dbetabinom(x = m, size = c, prob = vb[v_b], rho = betab_rho, log = TRUE)
		}
		LL = sum(apply(loglik, 1, function(x) { log(sum(exp(x))) }))
		if (LL>LL_b) {
			LL_b = LL
			params = list(LL = LL,
				      r_b = rb,
				      p_bj = p_bj)
		}
		p_bjr[[r]] = p_bj
	}
	params = list(LL = params$LL,
		      r_b = params$r_b,
		      p_bj = params$p_bj,
		      p_bjr = p_bjr)
	return(invisible(params))
}
