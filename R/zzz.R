.MMEnv <- new.env()

n_run = 100
betab_rho = 1e-4
max_iter = 1000
min_dist = 1e-4
vb = 2^(-2:-8)
vb2 <- function(nb) { rep(.MMEnv$vb[1:nb], each = 2) }

assign("n_run", n_run, envir=.MMEnv)
assign("betab_rho", betab_rho, envir=.MMEnv)
assign("max_iter", max_iter, envir=.MMEnv)
assign("min_dist", min_dist, envir=.MMEnv)
assign("vb", vb, envir=.MMEnv)
assign("vb2", vb2, envir=.MMEnv)

rm(list=c("n_run", "betab_rho", "max_iter", "min_dist", "vb", "vb2"))
