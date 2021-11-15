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
library('ggsignif')
library('diverse')

'prune_segments' <- function(x, n = 10)
{
	cnm = matrix(NA, nrow = nrow(x), ncol = nrow(x))
	for (j in 1:nrow(x)) {
		cnm[,j] = abs(2^x[j,"Log2Ratio"] - 2^x[,"Log2Ratio"])
	}
	cnt = hclust(as.dist(cnm), "average")
	cnc = cutree(tree = cnt, k = n)
	for (j in unique(cnc)) {
		indx = which(cnc==j)
		if (length(indx)>2) {
			mcl = mean(x[indx,"Log2Ratio"])
			scl = sd(x[indx,"Log2Ratio"])
			ind = which(x[indx,"Log2Ratio"]<(mcl+1.96*scl) & x[indx,"Log2Ratio"]>(mcl-1.96*scl))
			x[indx[ind],"Log2Ratio"] = mean(x[indx[ind],"Log2Ratio"])
		} else {
			x[indx,"Log2Ratio"] = mean(x[indx,"Log2Ratio"])
		}
	}
	return(invisible(x))
}

'absolute_copynumber' <- function(x, alpha, psi, c = 1)
{
	y = (1/alpha) * ((2^(x/c))*((psi*alpha) + 2*(1-alpha)) - 2*(1-alpha))
	return(invisible(y))
}

data("r(x)")
data("hg19")

data_ = data %>%
	dplyr::mutate(Position = round(.5*(Start + End))) %>%
	dplyr::select(Chromosome, Position, `01-T`) %>%
	dplyr::rename(Log2Ratio = `01-T`) %>%
	base::as.data.frame()
segmented = pcf(data = winsorize(data = data_, method = "mad", tau = 2.5, k = 25, verbose = FALSE), kmin = 50, gamma = 50, fast = FALSE, verbose = FALSE)[,2:7,drop = FALSE]
colnames(segmented) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
segmented = prune_segments(x = segmented, n = 10)
smoothed = winsorize(data = data_[,c("Chromosome","Position","Log2Ratio")], tau = 2.5, k = 15, verbose = FALSE)
colnames(smoothed) = c("Chromosome", "Position", "Log2Ratio")
end = NULL
for (j in 1:22) {
	end = c(end, max(hg19$End[hg19$Chromosome == j]))
}
end = cumsum(end)
start = rep(0, 22)
start[2:22] = end[1:21]+1
for (j in 1:22) {
	segmented[segmented[,"Chromosome"]==j,"Start"] = segmented[segmented[,"Chromosome"]==j,"Start"] + start[j]
	segmented[segmented[,"Chromosome"]==j,"End"] = segmented[segmented[,"Chromosome"]==j,"End"] + start[j]
	smoothed[smoothed[,"Chromosome"]==j,"Position"] = smoothed[smoothed[,"Chromosome"]==j,"Position"] + start[j]
}
col = rep("grey75", nrow(smoothed))
pdf(file = "01-T.pdf", width = 13, height=4)
par(mar=c(5, 5, 4, 2)+.1)
plot(x = smoothed[,"Position"], y = smoothed[,"Log2Ratio"], type = "p", pch = ".", cex = 1, col = col, axes = FALSE, frame = FALSE, xlab = "", ylab = "", main = "", ylim = c(-2.5,2))
for (j in 1:nrow(segmented)) {
	lines(x = c(segmented[j,"Start"], segmented[j,"End"]), y = rep(segmented[j,"Log2Ratio"], 2), lty = 1, lwd = 2.75, col = "red")
}
axis(side = 2, at = c(-2,-1,0,1,2), labels = c("-2.0","-1.0","0.0","1.0","2.0"), las = 1, cex.axis = 1, lwd = 1.5, lwd.ticks = 1.25, line = -.5)
axis(side = 1, at = c(start, end[length(end)]), labels = rep("", length(start)+1), tcl = .35, lwd = 1.5, lwd.ticks = 1.25)
axis(side = 1, at = .5*(start+end), labels = c(1:22), tcl = -.5, lwd = 0, lwd.ticks = 1.25)
mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
dev.off()

resi = matrix(NA, nrow = 100, ncol = 100)
alpha = seq(from = .1, to = 1, length = 100)
psi = seq(from = 1.5, to = 4, length = 100)
for (i in 1:length(alpha)) {
	for (j in 1:length(psi)) {
		abs_cn = absolute_copynumber(x = segmented$Log2Ratio, alpha = alpha[i], psi = psi[j], c = 1)
		resi[j,i] = sum(abs(round(abs_cn) - abs_cn) * (segmented$End - segmented$Start)/sum(segmented$End - segmented$Start))
	}
}

pdf(file = "01-T(2).pdf")
image(x = psi, y = alpha, z = resi,
      col = rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256)),
      las = 1,
      xlab = expression(psi),
      ylab = expression(alpha),
      main = "MOS01-T")
abline(h = 0.41, col = "green", lty = 2)
abline(v = 2.67, col = "green", lty = 2)
dev.off()

abs_cn = round(absolute_copynumber(x = segmented$Log2Ratio, alpha = .41, psi = 2.67, c = 1))
index = abs_cn!=2
prop_alt = segmented[index,,drop=FALSE]
sum(prop_alt$End - prop_alt$Start)/sum(segmented$End - segmented$Start)

data_ = data %>%
	dplyr::mutate(Position = round(.5*(Start + End))) %>%
	dplyr::select(Chromosome, Position, `08-T1`) %>%
	dplyr::rename(Log2Ratio = `08-T1`) %>%
	base::as.data.frame()
segmented = pcf(data = winsorize(data = data_, method = "mad", tau = 2.5, k = 25, verbose = FALSE), kmin = 50, gamma = 50, fast = FALSE, verbose = FALSE)[,2:7,drop = FALSE]
colnames(segmented) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
segmented = prune_segments(x = segmented, n = 17)
smoothed = winsorize(data = data_[,c("Chromosome","Position","Log2Ratio")], tau = 2.5, k = 15, verbose = FALSE)
colnames(smoothed) = c("Chromosome", "Position", "Log2Ratio")
end = NULL
for (j in 1:22) {
	end = c(end, max(hg19$End[hg19$Chromosome == j]))
}
end = cumsum(end)
start = rep(0, 22)
start[2:22] = end[1:21]+1
for (j in 1:22) {
	segmented[segmented[,"Chromosome"]==j,"Start"] = segmented[segmented[,"Chromosome"]==j,"Start"] + start[j]
	segmented[segmented[,"Chromosome"]==j,"End"] = segmented[segmented[,"Chromosome"]==j,"End"] + start[j]
	smoothed[smoothed[,"Chromosome"]==j,"Position"] = smoothed[smoothed[,"Chromosome"]==j,"Position"] + start[j]
}
col = rep("grey75", nrow(smoothed))
pdf(file = "08-T1.pdf", width = 13, height=4)
par(mar=c(5, 5, 4, 2)+.1)
plot(x = smoothed[,"Position"], y = smoothed[,"Log2Ratio"], type = "p", pch = ".", cex = 1, col = col, axes = FALSE, frame = FALSE, xlab = "", ylab = "", main = "", ylim = c(-2.5,2))
for (j in 1:nrow(segmented)) {
	lines(x = c(segmented[j,"Start"], segmented[j,"End"]), y = rep(segmented[j,"Log2Ratio"], 2), lty = 1, lwd = 2.75, col = "red")
}
axis(side = 2, at = c(-2,-1,0,1,2), labels = c("-2.0","-1.0","0.0","1.0","2.0"), las = 1, cex.axis = 1, lwd = 1.5, lwd.ticks = 1.25, line = -.5)
axis(side = 1, at = c(start, end[length(end)]), labels = rep("", length(start)+1), tcl = .35, lwd = 1.5, lwd.ticks = 1.25)
axis(side = 1, at = .5*(start+end), labels = c(1:22), tcl = -.5, lwd = 0, lwd.ticks = 1.25)
mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
dev.off()


resi = matrix(NA, nrow = 100, ncol = 100)
alpha = seq(from = .1, to = 1, length = 100)
psi = seq(from = 1.5, to = 4, length = 100)
for (i in 1:length(alpha)) {
	for (j in 1:length(psi)) {
		abs_cn = absolute_copynumber(x = segmented$Log2Ratio, alpha = alpha[i], psi = psi[j], c = 1)
		resi[j,i] = sum(abs(round(abs_cn) - abs_cn) * (segmented$End - segmented$Start)/sum(segmented$End - segmented$Start))
	}
}

pdf(file = "08-T1(2).pdf")
image(x = psi, y = alpha, z = resi,
      col = rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256)),
      las = 1,
      xlab = expression(psi),
      ylab = expression(alpha),
      main = "MOS08-T1")
abline(h = 0.36, col = "green", lty = 2)
abline(v = 2.15, col = "green", lty = 2)
dev.off()

abs_cn = round(absolute_copynumber(x = segmented$Log2Ratio, alpha = .36, psi = 2.15, c = 1))
index = abs_cn!=2
prop_alt = segmented[index,,drop=FALSE]
sum(prop_alt$End - prop_alt$Start)/sum(segmented$End - segmented$Start)


data_ = data %>%
	dplyr::mutate(Position = round(.5*(Start + End))) %>%
	dplyr::select(Chromosome, Position, `08-T2`) %>%
	dplyr::rename(Log2Ratio = `08-T2`) %>%
	base::as.data.frame()
segmented = pcf(data = winsorize(data = data_, method = "mad", tau = 2.5, k = 25, verbose = FALSE), kmin = 50, gamma = 50, fast = FALSE, verbose = FALSE)[,2:7,drop = FALSE]
colnames(segmented) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
segmented = prune_segments(x = segmented, n = 10)
smoothed = winsorize(data = data_[,c("Chromosome","Position","Log2Ratio")], tau = 2.5, k = 15, verbose = FALSE)
colnames(smoothed) = c("Chromosome", "Position", "Log2Ratio")
end = NULL
for (j in 1:22) {
	end = c(end, max(hg19$End[hg19$Chromosome == j]))
}
end = cumsum(end)
start = rep(0, 22)
start[2:22] = end[1:21]+1
for (j in 1:22) {
	segmented[segmented[,"Chromosome"]==j,"Start"] = segmented[segmented[,"Chromosome"]==j,"Start"] + start[j]
	segmented[segmented[,"Chromosome"]==j,"End"] = segmented[segmented[,"Chromosome"]==j,"End"] + start[j]
	smoothed[smoothed[,"Chromosome"]==j,"Position"] = smoothed[smoothed[,"Chromosome"]==j,"Position"] + start[j]
}
col = rep("grey75", nrow(smoothed))
pdf(file = "08-T2.pdf", width = 13, height=4)
par(mar=c(5, 5, 4, 2)+.1)
plot(x = smoothed[,"Position"], y = smoothed[,"Log2Ratio"], type = "p", pch = ".", cex = 1, col = col, axes = FALSE, frame = FALSE, xlab = "", ylab = "", main = "", ylim = c(-2.5,2))
for (j in 1:nrow(segmented)) {
	lines(x = c(segmented[j,"Start"], segmented[j,"End"]), y = rep(segmented[j,"Log2Ratio"], 2), lty = 1, lwd = 2.75, col = "red")
}
axis(side = 2, at = c(-2,-1,0,1,2), labels = c("-2.0","-1.0","0.0","1.0","2.0"), las = 1, cex.axis = 1, lwd = 1.5, lwd.ticks = 1.25, line = -.5)
axis(side = 1, at = c(start, end[length(end)]), labels = rep("", length(start)+1), tcl = .35, lwd = 1.5, lwd.ticks = 1.25)
axis(side = 1, at = .5*(start+end), labels = c(1:22), tcl = -.5, lwd = 0, lwd.ticks = 1.25)
mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
dev.off()

data_ = data %>%
	dplyr::mutate(Position = round(.5*(Start + End))) %>%
	dplyr::select(Chromosome, Position, `31-T`) %>%
	dplyr::rename(Log2Ratio = `31-T`) %>%
	base::as.data.frame()
segmented = pcf(data = winsorize(data = data_, method = "mad", tau = 2.5, k = 25, verbose = FALSE), kmin = 50, gamma = 50, fast = FALSE, verbose = FALSE)[,2:7,drop = FALSE]
colnames(segmented) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
segmented = prune_segments(x = segmented, n = 10)
smoothed = winsorize(data = data_[,c("Chromosome","Position","Log2Ratio")], tau = 2.5, k = 15, verbose = FALSE)
colnames(smoothed) = c("Chromosome", "Position", "Log2Ratio")
end = NULL
for (j in 1:22) {
	end = c(end, max(hg19$End[hg19$Chromosome == j]))
}
end = cumsum(end)
start = rep(0, 22)
start[2:22] = end[1:21]+1
for (j in 1:22) {
	segmented[segmented[,"Chromosome"]==j,"Start"] = segmented[segmented[,"Chromosome"]==j,"Start"] + start[j]
	segmented[segmented[,"Chromosome"]==j,"End"] = segmented[segmented[,"Chromosome"]==j,"End"] + start[j]
	smoothed[smoothed[,"Chromosome"]==j,"Position"] = smoothed[smoothed[,"Chromosome"]==j,"Position"] + start[j]
}
col = rep("grey75", nrow(smoothed))
pdf(file = "31-T.pdf", width = 13, height=4)
par(mar=c(5, 5, 4, 2)+.1)
plot(x = smoothed[,"Position"], y = smoothed[,"Log2Ratio"], type = "p", pch = ".", cex = 1, col = col, axes = FALSE, frame = FALSE, xlab = "", ylab = "", main = "", ylim = c(-2.5,2))
for (j in 1:nrow(segmented)) {
	lines(x = c(segmented[j,"Start"], segmented[j,"End"]), y = rep(segmented[j,"Log2Ratio"], 2), lty = 1, lwd = 2.75, col = "red")
}
axis(side = 2, at = c(-2,-1,0,1,2), labels = c("-2.0","-1.0","0.0","1.0","2.0"), las = 1, cex.axis = 1, lwd = 1.5, lwd.ticks = 1.25, line = -.5)
axis(side = 1, at = c(start, end[length(end)]), labels = rep("", length(start)+1), tcl = .35, lwd = 1.5, lwd.ticks = 1.25)
axis(side = 1, at = .5*(start+end), labels = c(1:22), tcl = -.5, lwd = 0, lwd.ticks = 1.25)
mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
dev.off()
