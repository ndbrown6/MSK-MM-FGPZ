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
library('Hmisc')
library('copynumber')

registerDoMC(8)

hex_cols = c("#C1272D",
	     "#377EB8",
	     "#01A99D",
	     "#F9ED7D",
	     "#F49C45")

'scientific_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^", scales::scientific_format()(x)), , fixed = TRUE))
}

data("vb=n")
m = data %>% .[["N_Alt"]]
c = data %>% .[["N_Total"]]

LL0 = MixtureBB(m = m, c = c, nb = 5)

# test asymmetric model with actual data vb = n and asymmetry a = 0.5->2
n = 250
ai = seq(from = 0.5, to = 2, length = n)
aj = seq(from = 0.5, to = 2, length = n)

LL = foreach(i=1:n) %dopar% {
	print(i)
	LL = vector(mode="numeric", length = n)
	for (j in 1:n) {
		# asymmetry of all branches
		a = rep(c(LL0$vb/.MMEnv$vb[1:5]), each = 2)
		a[seq(from = 2, to = 10, by = 2)] = 1
		
		# first cell division
		a[1] = ai[i]
		
		# second cell division
		a[3] = aj[j]
		
		LL[j] = AsymmLL(m = m, c = c, a = a, nb = 5)
	}
	return(invisible(LL))
}
save(list=ls(all=TRUE), file = "LL_Diff_1_2.RData")

LL = do.call(rbind, LL)

pdf(file = "LL_Diff_1_2.pdf", width = 6, height = 6)
par(mar=c(6.1, 6.5, 4.1, 3.1))
image(ai*.MMEnv$vb[1]*100, aj*.MMEnv$vb[2]*100, LL,
      xlab = expression("VAF of 1"^st~"cell division (%)"),
      ylab = expression("VAF of 2"^nd~"cell division (%)"),
      las=1, col = hcl.colors(55, "YlOrRd", rev = TRUE), useRaster = FALSE)
abline(v = .MMEnv$vb[1]*100, h = .MMEnv$vb[2]*100, lty = 2)
axis(side = 3, at = seq(from = 0.5, to = 2, length = 4)*.MMEnv$vb[1]*100, labels = seq(from = .5, to = 2, length = 4), cex.axis = .75, tcl = -.5)
axis(side = 3, at = seq(from = 0.75, to = 1.75, length = 3)*.MMEnv$vb[1]*100, labels = rep(" ", 3), cex.axis = .75, tcl = -.35)
axis(side = 4, at = seq(from = 0.5, to = 2, length = 4)*.MMEnv$vb[2]*100, labels = seq(from = .5, to = 2, length = 4), cex.axis = .75, tcl = -.5, las = 1)
axis(side = 4, at = seq(from = 0.75, to = 1.75, length = 3)*.MMEnv$vb[2]*100, labels = rep(" ", 3), cex.axis = .75, tcl = -.35)
contour(ai*.MMEnv$vb[1]*100, aj*.MMEnv$vb[2]*100, LL, add = TRUE, nlevels = 5, col = "grey10")
box()
dev.off()
