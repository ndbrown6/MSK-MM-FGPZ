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

registerDoMC(24)

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

save(list=ls(all=TRUE), file = "LL_Diff_4_5.RData")

LL = do.call(rbind, LL)
index = ai>=.5 & ai<=2
pdf(file = "LL_Diff_4_5.pdf", width = 5, height = 5.5)
par(mar=c(6.1, 6.5, 4.1, 1.1))
image(ai[index]*.MMEnv$vb[4]*100, aj[index]*.MMEnv$vb[5]*100, LL[index,index],
      xlab = expression("Expected VAF of 4"^th~"cell division"),
      ylab = expression("Expected VAF of 5"^th~"cell division"),
      las=1, col = hcl.colors(35, "YlOrRd", rev = TRUE))
abline(v = .MMEnv$vb[4]*100, h = .MMEnv$vb[5]*100, lty = 2)
contour(ai[index]*.MMEnv$vb[4]*100, aj[index]*.MMEnv$vb[5]*100, LL[index,index], add=TRUE, nlevels=15, col = "grey10")
dev.off()
