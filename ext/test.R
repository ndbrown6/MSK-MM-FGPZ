#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
library('mosaicm')

# test functions with simulated data at v_b = 1
data("vb=1")

ll = 0
for (i in 1:n) {
	ll = ll + log(TotalLL(m[i], c[i]))
}

# test asymmetric functions with simulated data at v_b = 1
data("vb=1")

ll = vector(mode = "numeric", length = 100)
a = seq(from = 1, to = 4, length = 100)
for (i in 1:100) {
	for (j in 1:n) {
		ll[i] = ll[i] + log(TotalLL(m[j], c[j], a = a[i], b = 2))
	}
}

# test asymmetric functions with simulated data at v_b = 0
data("vb=0")
m = data_$m
c = data_$c
n = nrow(data_)

'AsymmLL' <- function(m_j, c_j, v_b, a)
{
	ll = VGAM::dbetabinom(x = m_j, size = c_j, prob = v_b*a, rho = 1e-2)
	return(invisible(ll))
}

'TotalLL' <- function(m, c, a = NULL, b = NULL)
{
	a_b = rep(1, 15)
	a_b[b] = a
	ll = sum(AsymmLL(m, c, .MMEnv$v_b, a_b))
	return(invisible(ll))
}

ll = matrix(0, 50, 50)
a1 = seq(from = 1, to = 4, length = 50)
a2 = seq(from = 1, to = 8, length = 50)
for (i in 1:50) {
	for (j in 1:50) {
		ll0 = ll1 = 0
		for (ii in 1:n) {
			ll0 = ll0 + log(TotalLL(m[ii], c[ii], a = a1[i], b = 2))
			ll1 = ll1 + log(TotalLL(m[ii], c[ii], a = a2[j], b = 4))
		}
		ll[i,j] = ll0 - ll1
	}
}
