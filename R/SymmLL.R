'SymmLL' <- function(m, c)
{
	ll = AsymmLL(m = m, c = c, a = rep(1, .MMEnv$nb))
	return(invisible(ll))
}
