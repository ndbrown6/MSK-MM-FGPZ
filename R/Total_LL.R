'Total_LL' <- function(m, c, n, b, method = "binomial", log = TRUE)
{
	if (log) {
		ll = 0
		for (i in 1:n) {
			ll = ll + Symm_LL(m[i], c[i], b, method = method, log = log)
		}
	} else {
		ll = 1
		for (i in 1:n) {
			ll = ll * Symm_LL(m[i], c[i], b, method = method, log = log)
		}
	}
	return(invisible(ll))
}
