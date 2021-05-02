.MMEnv <- new.env()

v_b = c(.25, .125, .0625, .03125)
assign("v_b", v_b, envir=.MMEnv)
n_b = c(2, 4, 8, 16)
assign("n_b", n_b, envir=.MMEnv)

rm(v_b, n_b)
