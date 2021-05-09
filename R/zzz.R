.MMEnv <- new.env()

n_run = 100
betab_rho = 1e-3
min_vafdist = 0
max_iter = 100
min_dist = 0.0001
vb = list(c(0.25, 0.125, 0.0625, .03125), c(rep(.25, 1), rep(.125, 2), rep(.0625, 4), rep(.03125, 8)))
sb = list(rep(1, 4), rep(1, 15))

assign("n_run", n_run, envir=.MMEnv)
assign("betab_rho", betab_rho, envir=.MMEnv)
assign("min_vafdist", min_vafdist, envir=.MMEnv)
assign("max_iter", max_iter, envir=.MMEnv)
assign("min_dist", min_dist, envir=.MMEnv)
assign("vb", vb, envir=.MMEnv)
assign("sb", sb, envir=.MMEnv)

rm(list=c("n_run", "betab_rho", "min_vafdist", "max_iter", "min_dist", "vb", "sb"))
