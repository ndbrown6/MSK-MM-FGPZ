.GapEnv <- new.env()
vN = c("c_ind", "c_is", "c_if", "c_chr", "c_len",
	   "c_lrrBP", "c_bafBP", "c_lrr", "c_baf", "c_bafH",
	   "c_nC", "c_nA", "c_nS", "c_CN", "c_BA",
	   "c_conf", "c_sb", "c_CN1", "c_BA1", "c_rnk1")
for (i in 1:length(vN)) {
	assign(vN[i], i, envir=.GapEnv)
}