#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))

'fsq' <- function(alpha, beta, s_q, q_t)
{
	return(invisible(((alpha*s_q) + beta)/((alpha*q_t)+(1-alpha)*2)))
}

beta = 1 / seq(from = 2, to = 7, by = 1)
fsq_b = 2^seq(from = -2, to = -7, by = -1)
fsq_t = list()
alpha = seq(from = 0.2, to = .5, by = .1)
for (i in 1:length(alpha)) {
	fsq_t[[i]] = fsq(alpha = alpha[i], beta = beta, s_q = 1, q_t = 3)
}

pdf(file = "fsq_tumor_blood.pdf", width = 5, height = 5)
plot(0, 0, type = "n", xlab = expression(beta), ylab = "Tumor / Blood VAF Ratio", xlim = c(0.05, .5), ylim=c(0.5, 20), main = "", las = 1, log = "y")
for (i in 1:length(alpha)) {
	points(beta, fsq_t[[i]]/fsq_b, type = "l", col = "#772953")
}
abline(h = 1.5, col = "#B6AFA8", lwd = 2, lty = 1)
rect(xleft = .2, xright = .22, ybottom = 3, ytop = 20, col = "white", border = "white")
text(x = .21, y = 5, labels = "0.2", col = "#772953", cex = .55)
text(x = .21, y = 6.25, labels = "0.3", col = "#772953", cex = .55)
text(x = .21, y = 7.5, labels = "0.4", col = "#772953", cex = .55)
text(x = .21, y = 8.5, labels = "0.5", col = "#772953", cex = .55)
points(x = beta, y = rep(.5, length(beta)), pch = 25, col = "goldenrod3", bg = "goldenrod3", cex = .75)
box(lwd = 2)
dev.off()

beta_b = 1 / seq(from = 2, to = 7, by = 1)
fsq_b = 2^seq(from = -2, to = -7, by = -1)
fsq_t = list()
alpha = .2
beta_t = seq(from = .01, to = .5, by = .1)
for (i in 1:length(beta_t)) {
	fsq_t[[i]] = fsq(alpha = alpha, beta = beta_t[i], s_q = 1, q_t = 3)
}

pdf(file = "fsq_tumor_blood_noB.pdf", width = 5, height = 5)
plot(0, 0, type = "n", xlab = expression(beta^b), ylab = "Tumor / Blood VAF Ratio", xlim = c(0.05, .5), ylim=c(0.5, 20), main = "", las = 1, log = "y")
for (i in 1:length(beta_t)) {
	points(beta_b, fsq_t[[i]]/fsq_b, type = "l", col = "#772953")
}
abline(h = 1.5, col = "#B6AFA8", lwd = 2, lty = 1)
rect(xleft = .2, xright = .22, ybottom = 2.2, ytop = 20, col = "white", border = "white")
text(x = .21, y = 2.6, labels = "0.01", col = "#772953", cex = .55)
text(x = .21, y = 4.0, labels = "0.10", col = "#772953", cex = .55)
text(x = .21, y = 5.0, labels = "0.20", col = "#772953", cex = .55)
text(x = .21, y = 6.5, labels = "0.30", col = "#772953", cex = .55)
text(x = .215, y = 8.0, labels = "0.40", col = "#772953", cex = .55)
points(x = beta, y = rep(.5, length(beta)), pch = 25, col = "goldenrod3", bg = "goldenrod3", cex = .75)
box(lwd = 2)
dev.off()

