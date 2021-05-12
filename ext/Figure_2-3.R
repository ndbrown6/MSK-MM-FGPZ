#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
library('mosaicm')
library('ggplot2')
library('ggridges')
library('parallel')
library('foreach')
library('doMC')

registerDoMC(12)

hex_cols = c("#e41a1c",
	     "#377eb8",
	     "#4daf4a",
	     "#984ea3",
	     "#ff7f00")

'scientific_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^", scales::scientific_format()(x)), , fixed = TRUE))
}

data("vb=n")
m = data %>% .[["N_Alt"]]
c = data %>% .[["N_Total"]]

# test asymmetric model with actual data vb = n and asymmetry a = 0->3
n = 300
ai = seq(from = 0, to = 3, length = n)
aj = seq(from = 0, to = 3, length = n)

LL = foreach(i=1:n) %dopar% {
	LL = vector(mode="numeric", length = n)
	for (j in 1:n) {
		a_i = a_j = rep(1, 5)	
		a_i[4] = ai[i]
		a_j[5] = aj[j]
		LL[j] = AsymmLL(m = m, c = c, a = a_i, nb=5)$LL - AsymmLL(m = m, c = c, a = a_j, nb=5)$LL
	}
	return(invisible(LL))
}

save(list=ls(all=TRUE), file = "LL_Diff.RData")
