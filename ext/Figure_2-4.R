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

registerDoMC(4)

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

# test asymmetric model with actual data vb = n and asymmetry a = 0.5->2
n = 100
ai = seq(from = 0.5, to = 2, length = n)
aj = seq(from = 0.5, to = 2, length = n)

LL = foreach(i=1:n) %dopar% {
	print(i)
	LL = vector(mode="numeric", length = n)
	for (j in 1:n) {
		a = rep(1, 62)
		
		# third cell division
		a[7] = ai[i]
		
		# fourth cell division
		a[15] = aj[j]
		
		LL[j] = AsymmLL(m = m, c = c, a = a, nb=5)
	}
	return(invisible(LL))
}

#save(list=ls(all=TRUE), file = "LL_Diff_3_4.RData")

LL = do.call(rbind, LL)
index = ai>0 & ai<2
#pdf(file = "LL_Diff_3_4.pdf", width = 5, height = 5.5)
par(mar=c(6.1, 6.5, 4.1, 1.1))
image(ai[index]*.MMEnv$vb[3]*100, aj[index]*.MMEnv$vb[4]*100, LL[index,index],
      xlab = expression("Expected VAF of 3"^rd~"cell division"),
      ylab = expression("Expected VAF of 4"^th~"cell division"),
      las=1, col = hcl.colors(35, "YlOrRd", rev = TRUE))
abline(v = .MMEnv$vb[3]*100, h = .MMEnv$vb[4]*100, lty = 2)
contour(ai[index]*.MMEnv$vb[3]*100, aj[index]*.MMEnv$vb[4]*100, LL[index,index], add=TRUE, nlevels=5, col = "grey10")
#dev.off()
