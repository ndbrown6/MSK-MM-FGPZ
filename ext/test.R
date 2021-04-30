#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
library("mskmmfgpz")
data("vb=1")

lls = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll = vector(mode = "numeric", length = n)
	for (j in 1:n) {
		ll[i] = LL(m_j = m[j], c_j = c[j], b = i, alpha = .5)
	}
	lls[i] = sum(ll)
}

par(mar=c(6.1, 6.5, 4.1, 1.1))
barplot(height = lls, space = 0, names.arg = 1:4, las = 1, xlab = "")
mtext(side=1, text=expression(v[b]), line=4, cex=1.5)
mtext(side=2, text=expression(LL), line=4, cex=1.5)


alpha = seq(from = .5, to = 1, length = 100)
lls = vector(mode = "numeric", length = 4)
for (i in 1:4) {
	ll = vector(mode = "numeric", length = n)
	for (j in 1:n) {
		ll[i] = LL(m_j = m[j], c_j = c[j], b = i, alpha = .5)
	}
	ll_s[i] = sum(ll)
}

par(mar=c(6.1, 6.5, 4.1, 1.1))
barplot(height = ll_s, space = 0, names.arg = 1:4, las = 1, xlab = "")
mtext(side=1, text=expression(v[b]), line=4, cex=1.5)
mtext(side=2, text=expression(LL), line=4, cex=1.5)
