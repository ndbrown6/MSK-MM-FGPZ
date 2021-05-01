#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
library("mskmmfgpz")
data("vb=2")

'LL' <- function (m_j, c_j, b, rbsb) 
{
    ll = dbinom(x = m_j, size = c_j, prob = .MMEnv$v_b[b]) * rbsb *
         dbinom(x = c_j - m_j, size = c_j, prob = .MMEnv$v_b[b]) * (1-rbsb)
    return(invisible(ll))
}

'LL2' <- function (m_j, c_j, b, rbsb) 
{
    ll = dbeta(x=m_j/c_j, round(c_j*.MMEnv$v_b[b])+1, round(c_j*(1-.MMEnv$v_b[b]))+1) * rbsb *
         dbeta(x=(c_j-m_j)/c_j, round(c_j*.MMEnv$v_b[b])+1, round(c_j*(1-.MMEnv$v_b[b]))+1) * (1-rbsb)
    return(invisible(ll))
}



'ZZ' <- function(m_j, c_j, a1, a2) {
	ll = 0
	rbsb = c(a1, a2, .5, .5)
	for (b in 1:4) {
		ll = ll + LL(m_j, c_j, b, rbsb[b])
	}
	ll = sum(ll)/sum(rbsb)
}

a = b = seq(from = .01, to = .99, length = 100)
ab = matrix(NA, nrow = 100, ncol = 100)
for (ii in 1:100) {
	for (jj in 1:100) {
		ll = vector(mode="numeric", length = n)
		for (i in 1:n) {
			ll[i] = ZZ(m[i], c[i], a[ii], b[jj])
		}
		ll = sum(log(ll))
		ab[ii, jj] = ll
	}
}


par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(a, la, las = 1, xlab = "", ylab = "")
mtext(side=1, text=expression(alpha), line=4, cex=1.5)
mtext(side=2, text=expression(LL), line=4, cex=1.5)
