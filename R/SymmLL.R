'SymmLL' <- function(m, c, nb = 4)
{
	ll = AsymmLL(m = m, c = c, a = rep(1, nb), nb = nb)
	return(invisible(ll))
}
