b = 1:4
v_b = c(.25, .125, .0625, .03125)

'LL' <- function(m_j, c_j, n, b) {
	ll = 0
	for (j in 1:n) {
		ll = ll + dbinom(x = m_j[j], size = c_j[j], prob = v_b[b], log = TRUE) + dbinom(x = c_j[j] - m_j[j], size = c_j[j], prob = 1 - v_b[b], log = TRUE)
	}
	return(invisible(ll))
}

set.seed(0)
n = 50
c = round(runif(n = n, min = 100, max = 1000))
m = abs(round((c*v_b[4]) + runif(n = n, min = -10, max = 10)))

res = rep(NA, length = 4)
for (i in 1:4) {
	res[i] = LL(m, c, n, i)
}

par(mar=c(6.1, 6.5, 4.1, 1.1))
barplot(height = -res, space = 0, names.arg = 1:4, las = 1, xlab = "")
mtext(side=1, text=expression(v[b]), line=4, cex=1.5)
mtext(side=2, text=expression(LL), line=4, cex=1.5)
