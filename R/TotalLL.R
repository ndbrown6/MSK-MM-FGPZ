'TotalLL' <- function(m, c, a = NULL, b = NULL)
{
	if (is.null(a)) {
		ll = sum(SymmLL(m, c, .MMEnv$v_b))
	} else {
		a_b = rep(1, 15)
		a_b[b] = a
		ll = sum(AsymmLL(m, c, .MMEnv$v_b, a_b))
	}
	return(invisible(ll))
}
