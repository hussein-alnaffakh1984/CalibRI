# reflimR 1.1.0 point estimates and 95% CIs (rounding disabled), same call as the main worker.R.
# Usage: Rscript reflim_run.R <data_dir> <out_csv>
suppressMessages({ library(reflimR); library(data.table) })
a <- commandArgs(TRUE); DATA <- a[1]; OUT <- a[2]
meta <- fread(file.path(DATA, "meta.csv")); vals <- fread(cmd = paste("gzip -dc", shQuote(file.path(DATA, "values.csv.gz"))))
V <- split(vals$value, vals$ds_id)
res <- rbindlist(lapply(meta$ds_id, function(id) {
  x <- V[[as.character(id)]]
  q <- tryCatch({ z <- suppressWarnings(reflim(x, plot.it = FALSE, print.n = FALSE, apply.rounding = FALSE))
                  as.numeric(c(z$limits[c("lower.lim","upper.lim")], z$confidence.int[c("lower.lim.low","lower.lim.upp","upper.lim.low","upper.lim.upp")])) },
                error = function(e) rep(NA_real_, 6))
  data.table(ds_id = id, rl_LL = q[1], rl_UL = q[2], rl_CI95_LL_lo = q[3], rl_CI95_LL_hi = q[4], rl_CI95_UL_lo = q[5], rl_CI95_UL_hi = q[6])
}))
fwrite(res, OUT); cat("reflimR done:", nrow(res), "\n")
