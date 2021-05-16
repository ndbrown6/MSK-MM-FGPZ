.MMEnv <- new.env()

n_run = 100
betab_rho = 1e-4
max_iter = 1000
min_dist = 1e-4
vb = c(.25, .125, .0625, .03125, .015625, .0078125, .00390625)
vb15 <- function(nb) { rep(.MMEnv$vb[1:nb], times = 2^(1:nb)) }

assign("n_run", n_run, envir=.MMEnv)
assign("betab_rho", betab_rho, envir=.MMEnv)
assign("max_iter", max_iter, envir=.MMEnv)
assign("min_dist", min_dist, envir=.MMEnv)
assign("vb", vb, envir=.MMEnv)
assign("vb15", vb15, envir=.MMEnv)

rm(list=c("n_run", "betab_rho", "max_iter", "min_dist", "vb", "vb15"))
