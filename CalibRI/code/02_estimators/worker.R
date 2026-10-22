# Indirect-RI coverage audit worker.
# Usage: Rscript worker.R <mode: fast|boot> <data_dir> <out_dir> <shard> <n_shards> <cores> <hours_limit> [max_ds]
suppressMessages({ library(refineR); library(reflimR); library(parallel); library(data.table); library(future) })
plan(sequential)
a <- commandArgs(TRUE)
MODE <- a[1]; DATA <- a[2]; OUT <- a[3]; SHARD <- as.integer(a[4]); NSH <- as.integer(a[5])
CORES <- as.integer(a[6]); HOURS <- as.numeric(a[7]); MAXDS <- if(length(a) >= 8) as.integer(a[8]) else NA
T0 <- Sys.time(); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

meta <- fread(file.path(DATA, "meta.csv"))
vals <- fread(cmd = paste("gzip -dc", shQuote(file.path(DATA, "values.csv.gz"))))
V <- split(vals$value, vals$ds_id)

if (MODE == "fast") {
  todo <- meta$ds_id
} else {  # bootstrap subset: simulated TSH n=500 (50 per fraction) + NHANES TSH n=500 (20 per fraction)
  todo <- c(meta[source=="sim" & analyte=="TSH" & n_nominal==500 & rep<=50]$ds_id,
            meta[source=="nhanes" & analyte=="TSH" & n_nominal==500 & rep<=20]$ds_id)
}
todo <- sort(todo); todo <- todo[(seq_along(todo) - 1) %% NSH + 1 == SHARD]
if (!is.na(MAXDS)) todo <- head(todo, MAXDS)

resfile <- file.path(OUT, sprintf("res_%s_shard%02d.csv", MODE, SHARD))
done <- if (file.exists(resfile)) fread(resfile)$ds_id else integer(0)
todo <- setdiff(todo, done)
cat(sprintf("[%s] mode=%s shard=%d/%d: %d datasets to run (%d already done)\n", format(Sys.time()), MODE, SHARD, NSH, length(todo), length(done)))

na6 <- rep(NA_real_, 6)
run_one <- function(id) {
  x <- V[[as.character(id)]]; t1 <- Sys.time(); out <- list(ds_id = id)
  if (MODE == "fast") {
    r <- tryCatch({ f <- suppressMessages(findRI(x, NBootstrap = 0))
                    g90 <- getRI(f, UMprop = 0.90); g95 <- getRI(f, UMprop = 0.95)
                    c(g90$PointEst, g90$UMLow, g90$UMHigh, g95$UMLow, g95$UMHigh) }, error = function(e) rep(NA_real_, 10))
    out[c("rf_LL","rf_UL","rf_UM90_LL_lo","rf_UM90_UL_lo","rf_UM90_LL_hi","rf_UM90_UL_hi",
          "rf_UM95_LL_lo","rf_UM95_UL_lo","rf_UM95_LL_hi","rf_UM95_UL_hi")] <- as.list(r)
    t2 <- Sys.time()
    q <- tryCatch({ z <- suppressWarnings(reflim(x, plot.it = FALSE, print.n = FALSE, apply.rounding = FALSE))
                    c(z$limits[c("lower.lim","upper.lim")], z$confidence.int[c("lower.lim.low","lower.lim.upp","upper.lim.low","upper.lim.upp")]) },
                  error = function(e) na6)
    out[c("rl_LL","rl_UL","rl_CI95_LL_lo","rl_CI95_LL_hi","rl_CI95_UL_lo","rl_CI95_UL_hi")] <- as.list(as.numeric(q))
    np <- quantile(x, c(.025, .975), names = FALSE); out$np_LL <- np[1]; out$np_UL <- np[2]
    out$sec_refineR <- as.numeric(difftime(t2, t1, units = "secs")); out$sec_total <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  } else {
    r <- tryCatch({ f <- suppressMessages(findRI(x, NBootstrap = 50, seed = id))
                    g <- getRI(f, CIprop = 0.95); g9 <- getRI(f, CIprop = 0.90)
                    c(g$PointEst, g$CILow, g$CIHigh, g9$CILow, g9$CIHigh) }, error = function(e) rep(NA_real_, 10))
    out[c("rfb_LL","rfb_UL","rfb_CI95_LL_lo","rfb_CI95_UL_lo","rfb_CI95_LL_hi","rfb_CI95_UL_hi",
          "rfb_CI90_LL_lo","rfb_CI90_UL_lo","rfb_CI90_LL_hi","rfb_CI90_UL_hi")] <- as.list(r)
    out$sec_total <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  }
  as.data.table(out)
}

# process in chunks so results are saved continuously and the time limit is respected
chunk <- CORES * 4
for (s in seq(1, length(todo), by = chunk)) {
  el <- as.numeric(difftime(Sys.time(), T0, units = "hours"))
  if (el > HOURS) { cat(sprintf("Time limit reached (%.2f h). Stopping; re-run the same command to resume.\n", el)); break }
  ids <- todo[s:min(s + chunk - 1, length(todo))]
  res <- rbindlist(mclapply(ids, run_one, mc.cores = CORES), fill = TRUE)
  fwrite(res, resfile, append = file.exists(resfile))
  n_done <- s + length(ids) - 1
  cat(sprintf("[%s] %d/%d done (%.2f h elapsed)\n", format(Sys.time(), "%H:%M:%S"), n_done, length(todo), as.numeric(difftime(Sys.time(), T0, units = "hours"))))
}
cat("Finished shard.\n")
