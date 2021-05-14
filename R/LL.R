'LL' <- function(m, c, nb)
{
	n_run = .MMEnv$n_run
	betab_rho = .MMEnv$betab_rho
	max_iter = .MMEnv$max_iter
	min_dist = .MMEnv$min_dist
	vb = .MMEnv$vb[1:nb]
	
	LL_b = -Inf
	p_bjr = list()
	for (r in 1:n_run) {
		vb_old = rep(-1, times = nb)
		vb_est = vb
		iter = 0
		while (iter<max_iter & sum(abs(vb_est-vb_old))>min_dist) {
			iter = iter + 1
			vb_old = vb_est
			p_bj = matrix(NA, nrow = length(m), ncol = nb)
			for (v_b in 1:nb) {
				p_bj[,v_b] = VGAM::dbetabinom(x = m, size = c, prob = vb_est[v_b], rho = betab_rho) * 2^v_b
			}
			p_bj = pmax(p_bj, 1e-20)
			p_bj = p_bj/apply(p_bj, 1, sum)
			for (v_b in 1:nb) {
				vb_est[v_b] = sum(m*p_bj[,v_b])/sum(c*p_bj[,v_b])
			}
		}
		w_b = apply(p_bj, 2, sum)/sum(p_bj)
		loglik = matrix(NA, nrow = length(m), ncol = nb)
		for (v_b in 1:nb) {
			loglik[,v_b] = log(w_b[v_b]) + log(2^v_b) + VGAM::dbetabinom(x = m, size = c, prob = vb_est[v_b], rho = betab_rho, log = TRUE)
		}
		LL = sum(apply(loglik, 1, function(x) { log(sum(exp(x))) }))
		if (LL>LL_b) {
			LL_b = LL
			params = list(LL = LL,
				      vb = vb_est,
				      wb = w_b,
				      p_bj = p_bj)
		}
		p_bjr[[r]] = p_bj
	}
	params = list(LL = params$LL,
		      v_b = params$vb,
		      w_b = params$wb,
		      p_bj = params$p_bj,
		      p_bjr = p_bjr)
	return(invisible(params))

}
