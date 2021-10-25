#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))

'fsq' <- function(alpha, beta, s_q, q_t)
{
	return(invisible(((alpha*s_q) + beta)/((alpha*q_t)+(1-alpha)*2)))
}


pdf(file = "fsq_by_alpha.pdf", width = 7, height = 7)
par(mfrow = c(2, 2))

alpha = seq(from = 0, to = .5, by = .01)
beta = 0.5
s_q = c(1, 1, 2, 1, 2, 3)
q_t = c(1, 2, 2, 3, 3, 3)
col = c("steelblue", "goldenrod3", "salmon")
lty = c(1, 2, 3)
plot(0, 0, type = "n", xlab = expression(alpha), ylab = expression("f"[sq]), xlim = c(0,.5), ylim=c(0,1), main = expression(beta~"= 0.5"), las = 1)
for (i in 1:length(s_q)) {
	f_sq = fsq(alpha = alpha, beta = beta, s_q = s_q[i], q_t = q_t[i])
	points(alpha, f_sq, type = "l", lwd = 3, col = col[q_t[i]], lty = lty[s_q[i]])
}
abline(v = .2, col = "#B6AFA8", lwd = 2, lty = 1)
abline(h = fsq(alpha = .2, beta = beta, s_q = 1, q_t = 3), col = "#B6AFA8", lwd = 2, lty = 1)
legend(x = 0, y = 1.05, legend = c(1, 2, 3), lty = 1, lwd = 2, col = col, box.lwd = NA, text.col = "grey30")
legend(x = 0, y = .75, legend = c(1, 2, 3), lty = lty, lwd = 2, box.lwd = NA, col = "grey30", text.col = "grey30")
box(lwd = 2)

alpha = seq(from = 0, to = .5, by = .01)
beta = 0.25
s_q = c(1, 1, 2, 1, 2, 3)
q_t = c(1, 2, 2, 3, 3, 3)
col = c("steelblue", "goldenrod3", "salmon")
lty = c(1, 2, 3)
plot(0, 0, type = "n", xlab = expression(alpha), ylab = expression("f"[sq]), xlim = c(0,.5), ylim=c(0,1), main = expression(beta~"= 0.25"), las = 1)
for (i in 1:length(s_q)) {
	f_sq = fsq(alpha = alpha, beta = beta, s_q = s_q[i], q_t = q_t[i])
	points(alpha, f_sq, type = "l", lwd = 3, col = col[q_t[i]], lty = lty[s_q[i]])
}
abline(v = .2, col = "#B6AFA8", lwd = 2, lty = 1)
abline(h = fsq(alpha = .2, beta = beta, s_q = 1, q_t = 3), col = "#B6AFA8", lwd = 2, lty = 1)
box(lwd = 2)

alpha = seq(from = 0, to = .5, by = .01)
beta = 0.125
s_q = c(1, 1, 2, 1, 2, 3)
q_t = c(1, 2, 2, 3, 3, 3)
col = c("steelblue", "goldenrod3", "salmon")
lty = c(1, 2, 3)
plot(0, 0, type = "n", xlab = expression(alpha), ylab = expression("f"[sq]), xlim = c(0,.5), ylim=c(0,1), main = expression(beta~"= 0.125"), las = 1)
for (i in 1:length(s_q)) {
	f_sq = fsq(alpha = alpha, beta = beta, s_q = s_q[i], q_t = q_t[i])
	points(alpha, f_sq, type = "l", lwd = 3, col = col[q_t[i]], lty = lty[s_q[i]])
}
abline(v = .2, col = "#B6AFA8", lwd = 2, lty = 1)
abline(h = fsq(alpha = .2, beta = beta, s_q = 1, q_t = 3), col = "#B6AFA8", lwd = 2, lty = 1)
box(lwd = 2)

plot(0, 0, type = "n", xlab = expression(alpha), ylab = expression("f"[sq]), xlim = c(0, .5), ylim=c(0, .5), main = "", las = 1)
alpha = seq(from = 0, to = .5, by = .01)
beta = seq(from = 0.01, to = .6, by = .1)
s_q = 1
q_t = 3
for (i in 1:length(alpha)) {
	for (j in 1:length(beta)) {
		f_sq = fsq(alpha = alpha, beta = beta[j], s_q = s_q, q_t = q_t)
		points(alpha, f_sq, type = "l", col = "#772953")
	}
}
abline(v = .2, col = "#B6AFA8", lwd = 2, lty = 1)
abline(h = fsq(alpha = .2, beta = .01, s_q = s_q, q_t = q_t), col = "#B6AFA8", lwd = 2, lty = 2)
rect(xleft = .35, xright = .4, ybottom = .1, ytop = .5, col = "white", border = "white")
text(x = .3775, y = .16, labels = "0.01", col = "#772953", cex = .75)
text(x = .3775, y = .2, labels = "0.10", col = "#772953", cex = .75)
text(x = .3775, y = .25, labels = "0.20", col = "#772953", cex = .75)
text(x = .3775, y = .29, labels = "0.30", col = "#772953", cex = .75)
text(x = .3775, y = .335, labels = "0.40", col = "#772953", cex = .75)
text(x = .3775, y = .375, labels = "0.50", col = "#772953", cex = .75)
box(lwd = 2)
dev.off()
