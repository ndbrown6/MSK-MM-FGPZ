.MMEnv <- new.env()

v_b = c(rep(.25, 1), rep(.125, 2), rep(.0625, 4), rep(.03125, 8))
assign("v_b", v_b, envir=.MMEnv)

rm(v_b)
