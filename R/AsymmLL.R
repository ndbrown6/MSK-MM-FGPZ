'AsymmLL' <- function(m, c, a, nb = 4)
{
	n_run = .MMEnv$n_run
	betab_rho = .MMEnv$betab_rho
	min_vafdist = .MMEnv$min_vafdist
	max_iter = .MMEnv$max_iter
	min_dist = .MMEnv$min_dist
	if (nb == 4) {
		vb = .MMEnv$vb[[1]]*a
		sb = .MMEnv$sb[[1]]
	} else if (nb == 15) {
		vb = .MMEnv$vb[[2]]*a
		sb = .MMEnv$sb[[2]]
	}

	LL_b = -Inf
	P_bj = list()
	for (r in 1:n_run) {
		rb_prev = rep(-1, times = nb)
		rb_est = rep(runif(n = 1, min = 1e-3, max = 1), times = nb)
		iter = 0
		while (iter<max_iter & sum(abs(rb_est-rb_prev))>min_dist) {
			iter = iter + 1
			rb_prev = rb_est
			p_bj = matrix(NA, nrow = length(m), ncol = nb)
			for (v_b in 1:nb) {
				if (nb == 4) {
					p_bj[,v_b] = VGAM::dbetabinom(x = m, size = c, prob = vb[v_b], rho = betab_rho)*rb_est[v_b]*sb[v_b]*2^v_b
				} else if (nb == 15) {
					p_bj[,v_b] = VGAM::dbetabinom(x = m, size = c, prob = vb[v_b], rho = betab_rho)*rb_est[v_b]*sb[v_b]
				}
			}
			p_bj = pmax(p_bj, 1e-20)
			p_bj = p_bj/apply(p_bj, 1, sum)
			Nb = rep(NA, nb)
			for (v_b in 1:nb) {
				Nb[v_b] = sum(m*p_bj[,v_b])/sum(c*p_bj[,v_b])
			}
			rb_est = Nb/sum(Nb)
		}
		wb = apply(p_bj, 2, sum)/sum(p_bj)
		loglik = matrix(NA, nrow = length(m), ncol = nb)
		for (v_b in 1:nb) {
			if (nb == 4) {
				loglik[,v_b] = log(wb[v_b]) + VGAM::dbetabinom(x = m, size = c, prob = vb[v_b], rho = betab_rho, log = TRUE) + log(rb_est[v_b]) + log(sb[v_b]) + 2^v_b
			} else {
				loglik[,v_b] = log(wb[v_b]) + VGAM::dbetabinom(x = m, size = c, prob = vb[v_b], rho = betab_rho, log = TRUE) + log(rb_est[v_b]) + log(sb[v_b])
			}
		}
		#LL = sum(apply(loglik, 1, function(x) { log(sum(exp(x - max(x)))) + max(x) }))
		LL = sum(apply(loglik, 1, function(x) { sum(exp(x))/sum(rb_est*sb) }))
		if (LL>LL_b) {
			LL_b = LL
			params = list(LL = LL,
				      r_b = rb_est,
				      w_b = wb,
				      p_bj = p_bj)
		}
		P_bj[[r]] = p_bj
	}
	params = list(LL = params$LL,
		      r_b = params$r_b,
		      w_b = params$w_b,
		      p_bj = params$p_bj,
		      p_bjr = P_bj)
	return(invisible(params))
}
