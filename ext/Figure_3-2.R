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

hex_cols = c("#C1272D",
	     "#377EB8",
	     "#01A99D",
	     "#F9ED7D",
	     "#F49C45")

n = 15

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

data("r(y)")

S = seq(from = 1, to = 20, by = 1)
nst = matrix(0, nrow = ncol(data)-2, ncol = length(S))
for (i in 1:(ncol(data)-2)) {
	segmented = pcf(data = winsorize(data = data %>%
					 	dplyr::select(all_of(c(1,2,i+2))) %>%
					 	base::data.frame(),
					 method = "mad", tau = 2.5, k = 25, verbose = FALSE), kmin = 50, gamma = 50, fast = FALSE, verbose = FALSE)[,2:7,drop = FALSE]
	colnames(segmented) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	segmented = prune_segments(x = segmented, n = n-5)

	for (ii in 1:length(S)) {
		for (jj in 1:22) {
			log2 = subset(segmented, Chromosome == jj)
			if (nrow(log2)>1) {
				for (kk in 1:(nrow(log2)-1)) {
					sx = log2[kk,"End"] - log2[kk,"Start"]
					sy = log2[kk+1,"End"] - log2[kk+1,"Start"]
					if (sx >= (S[ii]*10^6) & sy >= (S[ii]*10^6)) {
						if (log2[kk,"Log2Ratio"]!=log2[kk+1,"Log2Ratio"]) {
							nst[i,ii] = nst[i,ii] + 1
						}
					}
				}
			}
		}
	}
}
nst = t(nst)
colnames(nst) = colnames(data)[c(-1,-2)]
nst = dplyr::as_tibble(nst) %>%
      reshape2::melt() %>%
      dplyr::as_tibble() %>%
      dplyr::bind_cols(S = rep(S, times = ncol(data)-2)) %>%
      dplyr::rename(sample_name = variable,
		    lst = value)

data("r(x)")

data_ = data %>%
	dplyr::mutate(Position = round(.5*(Start + End))) %>%
	dplyr::select(Chromosome, Position, `01-T`) %>%
	dplyr::rename(Log2Ratio = `01-T`) %>%
	base::as.data.frame()
segmented = pcf(data = winsorize(data = data_, method = "mad", tau = 2.5, k = 25, verbose = FALSE), kmin = 50, gamma = 50, fast = FALSE, verbose = FALSE)[,2:7,drop = FALSE]
colnames(segmented) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
segmented = prune_segments(x = segmented, n = n)

S = seq(from = 1, to = 20, by = 1)
n_st = vector(mode = "numeric", length = length(S))
for (i in 1:length(S)) {
	for (j in 1:22) {
		log2 = subset(segmented, Chromosome == j)
		if (nrow(log2)>1) {
			for (k in 1:(nrow(log2)-1)) {
				sx = log2[k,"End"] - log2[k,"Start"]
				sy = log2[k+1,"End"] - log2[k+1,"Start"]
				if (sx >= (S[i]*10^6) & sy >= (S[i]*10^6)) {
					if (log2[k,"Log2Ratio"]!=log2[k+1,"Log2Ratio"]) {
						n_st[i] = n_st[i] + 1
					}
				}
			}
		}
	}
}

nst = nst %>%
      dplyr::bind_rows(
	      dplyr::tibble(sample_name = rep("01-T", length(n_st)),
			    lst = n_st,
			    S = S))

data_ = data %>%
	dplyr::mutate(Position = round(.5*(Start + End))) %>%
	dplyr::select(Chromosome, Position, `08-T1`) %>%
	dplyr::rename(Log2Ratio = `08-T1`) %>%
	base::as.data.frame()
segmented = pcf(data = winsorize(data = data_, method = "mad", tau = 2.5, k = 25, verbose = FALSE), kmin = 50, gamma = 50, fast = FALSE, verbose = FALSE)[,2:7,drop = FALSE]
colnames(segmented) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
segmented = prune_segments(x = segmented, n = n)

S = seq(from = 1, to = 20, by = 1)
n_st = vector(mode = "numeric", length = length(S))
for (i in 1:length(S)) {
	for (j in 1:22) {
		log2 = subset(segmented, Chromosome == j)
		if (nrow(log2)>1) {
			for (k in 1:(nrow(log2)-1)) {
				sx = log2[k,"End"] - log2[k,"Start"]
				sy = log2[k+1,"End"] - log2[k+1,"Start"]
				if (sx >= (S[i]*10^6) & sy >= (S[i]*10^6)) {
					if (log2[k,"Log2Ratio"]!=log2[k+1,"Log2Ratio"]) {
						n_st[i] = n_st[i] + 1
					}
				}
			}
		}
	}
}

nst = nst %>%
      dplyr::bind_rows(
	      dplyr::tibble(sample_name = rep("08-T1", length(n_st)),
			    lst = n_st,
			    S = S))

data_ = data %>%
	dplyr::mutate(Position = round(.5*(Start + End))) %>%
	dplyr::select(Chromosome, Position, `08-T2`) %>%
	dplyr::rename(Log2Ratio = `08-T2`) %>%
	base::as.data.frame()
segmented = pcf(data = winsorize(data = data_, method = "mad", tau = 2.5, k = 25, verbose = FALSE), kmin = 50, gamma = 50, fast = FALSE, verbose = FALSE)[,2:7,drop = FALSE]
colnames(segmented) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
segmented = prune_segments(x = segmented, n = n)

S = seq(from = 1, to = 20, by = 1)
n_st = vector(mode = "numeric", length = length(S))
for (i in 1:length(S)) {
	for (j in 1:22) {
		log2 = subset(segmented, Chromosome == j)
		if (nrow(log2)>1) {
			for (k in 1:(nrow(log2)-1)) {
				sx = log2[k,"End"] - log2[k,"Start"]
				sy = log2[k+1,"End"] - log2[k+1,"Start"]
				if (sx >= (S[i]*10^6) & sy >= (S[i]*10^6)) {
					if (log2[k,"Log2Ratio"]!=log2[k+1,"Log2Ratio"]) {
						n_st[i] = n_st[i] + 1
					}
				}
			}
		}
	}
}

nst = nst %>%
      dplyr::bind_rows(
	      dplyr::tibble(sample_name = rep("08-T2", length(n_st)),
			    lst = n_st,
			    S = S))

data_ = data %>%
	dplyr::mutate(Position = round(.5*(Start + End))) %>%
	dplyr::select(Chromosome, Position, `31-T`) %>%
	dplyr::rename(Log2Ratio = `31-T`) %>%
	base::as.data.frame()
segmented = pcf(data = winsorize(data = data_, method = "mad", tau = 2.5, k = 25, verbose = FALSE), kmin = 50, gamma = 50, fast = FALSE, verbose = FALSE)[,2:7,drop = FALSE]
colnames(segmented) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
segmented = prune_segments(x = segmented, n = n)

S = seq(from = 1, to = 20, by = 1)
n_st = vector(mode = "numeric", length = length(S))
for (i in 1:length(S)) {
	for (j in 1:22) {
		log2 = subset(segmented, Chromosome == j)
		if (nrow(log2)>1) {
			for (k in 1:(nrow(log2)-1)) {
				sx = log2[k,"End"] - log2[k,"Start"]
				sy = log2[k+1,"End"] - log2[k+1,"Start"]
				if (sx >= (S[i]*10^6) & sy >= (S[i]*10^6)) {
					if (log2[k,"Log2Ratio"]!=log2[k+1,"Log2Ratio"]) {
						n_st[i] = n_st[i] + 1
					}
				}
			}
		}
	}
}

nst = nst %>%
      dplyr::bind_rows(
	      dplyr::tibble(sample_name = rep("31-T", length(n_st)),
			    lst = n_st,
			    S = S))

plot_ = nst %>%
	dplyr::mutate(hrd = case_when(
		sample_name %in% c("01-T", "08-T1", "08-T2", "31-T") ~ "Yes",
		TRUE ~ "No"
	)) %>%
	ggplot(aes(x = S, y = lst, group = sample_name, color = hrd)) +
	geom_step(stat = "identity", direction = "vh", size = 1, alpha = .75) +
	geom_vline(xintercept = c(6,11), linetype = 3, size = .5, color = "goldenrod3") +
	scale_color_manual(values = hex_cols[2:1]) +
	xlab("\n\nSegment size (Mb)\n") +
	ylab("\nNumber of state transitions\n\n") +
	scale_y_log10(limits = c(10,75)) +
	annotation_logticks(side = "l") +
	theme_classic() +
	guides(color = guide_legend(title = "HRD"))

pdf(file = "nst.pdf", width = 6.5, height = 5)
print(plot_)
dev.off()

