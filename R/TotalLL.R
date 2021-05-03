'TotalLL' <- function(m, c, b, a = NULL, method = "binomial", log = TRUE)
{
	if (is.null(a)) {
		if (log) {
			ll = sum(SymmLL(m, c, b, method = method, log = TRUE))
		} else {
			ll = exp(sum(log(SymmLL(m, c, b, method = method, log = FALSE))))
		}
	} else {
		if (b==1) {
			if (a < 0 | a > .MMEnv$n_b[1]) {
				stop("a is outside of expected range!\n")
			}
		} else if (b==2) {
			if (a < 0 | a > .MMEnv$n_b[2]) {
				stop("a is outside of expected range!\n")
			}
		} else if (b==3) {
			if (a < 0 | a > .MMEnv$n_b[3]) {
				stop("a is outside of expected range!\n")
			}
		} else if (b==4) {
			if (a < 0 | a > .MMEnv$n_b[4]) {
				stop("a is outside of expected range!\n")
			}
		}
		if (log) {
			ll = sum(AsymmLL(m, c, b, a, method = method, log = TRUE))
		} else {
			ll = exp(sum(log(AsymmLL(m, c, b, a, method = method, log = FALSE))))
		}
	}
	return(invisible(ll))
}
