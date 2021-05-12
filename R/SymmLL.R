'SymmLL' <- function(m, c, nb)
{
	params = AsymmLL(m, c, a = rep(1, nb), nb = nb)
	return(invisible(params))
}
