.MMEnv <- new.env()

b = 1:4
v_b = c(.25, .125, .0625, .03125)
assign("b", b, envir=.MMEnv)
assign("v_b", v_b, envir=.MMEnv)
