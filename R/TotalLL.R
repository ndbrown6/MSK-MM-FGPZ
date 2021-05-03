'TotalLL' <- function(m, c, b, method = "binomial", log = TRUE)
{
	if (log) {
		ll = sum(SymmLL(m[i], c[i], b, method = method, log = TRUE))
	} else {
		ll = exp(sum(log(SymmLL(m[i], c[i], b, method = method, log = FALSE))))
	}
	return(invisible(ll))
}
