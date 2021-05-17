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

# test vb free model with actual data vb = n and nb = 1|2|3|4|5|6|7 cell generations
LL = vector(mode = "numeric", length = 7)
for (i in 1:7) {
	LL[i] = LL(m = m, c = c, nb = i)$LL
}

P = 1 - pchisq(2*(LL[2]-LL[1]),1)
P = 1 - pchisq(2*(LL[3]-LL[2]),1)
P = 1 - pchisq(2*(LL[4]-LL[3]),1)
P = 1 - pchisq(2*(LL[5]-LL[4]),1)
P = 1 - pchisq(2*(LL[6]-LL[5]),1)
P = 1 - pchisq(2*(LL[7]-LL[6]),1)
