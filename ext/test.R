#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
library("mskmmfgpz")
data("vb=1")

m = m + 100
m[m>c] = c[m>c]

alpha = seq(from = 1, to = 2, length = 100)
LP_A = vector(mode = "numeric", length = 100)
for (a in 1:length(alpha)) {
	ll = vector(mode = "numeric", length = n)
	for (i in 1:n) {
		ll[i] = LL(m_j = m[i], c_j = c[i], b = 1, alpha = alpha[a])
	}
	LP_A[a] = sum(log(ll))
}

par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(alpha, LP_A, las = 1, xlab = "", ylab = "")
mtext(side=1, text=expression(alpha), line=4, cex=1.5)
mtext(side=2, text=expression(LL), line=4, cex=1.5)
