#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
library("viridis")
library("pwr")

########## VAF versus coverage
n_len = 1000
p1 = seq(from = 0.01125, to = 0.25, length = n_len)
p2 = 10E-4
power = c(0.8, 0.9, 0.95, 0.99, 0.999)
n  = list()

for (i in 1:length(power)) {
	
	n[[i]] = vector(mode = "numeric", length = n_len)
	
	for (j in 1:n_len) {
		
		n[[i]][j] = pwr.p.test(h = ES.h(p1 = p1[j], p2 = p2), sig.level = 0.05, power = power[i], alternative = "greater")$n
	}
	
}

pdf(file = "Coverage_by_VAF_for_Power.pdf", width = 5, height = 5)
par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(1, 1, type = "n", xlim = c(0, .25)*100, ylim = c(1, 1000), xlab = "", ylab = "", las = 1, log = "y")
for (i in 1:length(power)) {
	points(p1*100, n[[i]], type = "l", col = viridis(n = length(power))[i], lwd = 3)
}
mtext(side = 1, text = "VAF (%)", line = 4, cex = 1)
mtext(side = 2, text = "Coverage", line = 5, cex = 1)
box(lwd = 2)
dev.off()

########## VAF versus power
n_len = 1000
p1 = seq(from = 0.001, to = 0.25, length = n_len)
p2 = 10E-4
n = c(20, 30, 40, 50, 100, 500, 1000)
power  = list()

for (i in 1:length(n)) {
	
	power[[i]] = vector(mode = "numeric", length = n_len)
	
	for (j in 1:n_len) {
		
		power[[i]][j] = pwr.p.test(h = ES.h(p1 = p1[j], p2 = p2), sig.level = 0.05, n = n[i], alternative = "greater")$power
	}
	
}

pdf(file = "Power_by_VAF_for_N.pdf", width = 5, height = 5)
par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type = "n", xlim = c(0, .25)*100, ylim = c(0, 1), xlab = "", ylab = "", las = 1)
for (i in 1:length(n)) {
	points(p1*100, power[[i]], type = "l", col = viridis(n = length(power))[i], lwd = 3)
}
mtext(side = 1, text = "VAF (%)", line = 4, cex = 1)
mtext(side = 2, text = "Probability", line = 5, cex = 1)
box(lwd = 2)
dev.off()

########## coverage versus power
n_len = 1000
p1 = 0.015
p2 = 10E-4
power = .8

n = pwr.p.test(h = ES.h(p1 = p1, p2 = p2), sig.level = 0.05, power = power, alternative = "greater")$n
power = pwr.p.test(h = ES.h(p1 = p1, p2 = p2), sig.level = 0.05, n = 250, alternative = "greater")$power