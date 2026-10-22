# ------------------------------------------------------------------------------------------
# Check: official rinet R package version 0.1.1 (the version used in CalibRI) versus the
# Python re-implementation used in the study.
# The current CRAN version of rinet loads its model with torch; version 0.1.1 loads the
# Keras model through reticulate. This script therefore installs 0.1.1 from the CRAN archive.
#
# Google Colab, R runtime. IMPORTANT: start from a fresh session
# (Runtime > Restart session, or a new notebook with Runtime type = R), then run this ONE cell.
# ------------------------------------------------------------------------------------------
options(repos = "https://cloud.r-project.org")
if ("rinet" %in% loadedNamespaces())
  stop("A previous rinet version is still loaded in memory. Use Runtime > Restart session, then run this cell again.")
if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
library(reticulate)
# Declare the Python packages BEFORE Python starts (reticulate creates an environment with them)
py_require(c("tensorflow", "keras", "scikit-learn", "scipy", "numpy"))

# Install rinet 0.1.1 from the CRAN archive (replaces any newer torch-based version)
has_rinet <- "rinet" %in% rownames(installed.packages())       # does not load the package
if (!has_rinet || utils::packageDescription("rinet")$Version != "0.1.1") {
  if (has_rinet) remove.packages("rinet")
  install.packages("https://cran.r-project.org/src/contrib/Archive/rinet/rinet_0.1.1.tar.gz",
                   repos = NULL, type = "source")
}
library(rinet)
stopifnot(!("model_version" %in% names(formals(predict_rinet))))   # must be the 0.1.1 code
cat("rinet version:", as.character(packageVersion("rinet")), "\n")
cat("Python used:", py_config()$python, "\n")
cat("tensorflow", py_eval("__import__('tensorflow').__version__"),
    "| keras", py_eval("__import__('keras').__version__"), "\n")

# 60 test datasets: log-normal healthy component + log-normal pathological component
set.seed(20261022)
make_ds <- function() {
  n <- sample(c(300, 500, 1000), 1); f <- sample(c(0.1, 0.3), 1); k <- round(n * f)
  sample(c(rlnorm(n - k, log(runif(1, 1, 3)), runif(1, 0.3, 0.6)),
           rlnorm(k, log(runif(1, 5, 15)), runif(1, 0.4, 0.8))))
}
datasets <- replicate(60, make_ds(), simplify = FALSE)

# (1) Official R package 0.1.1
r_res <- t(sapply(datasets, function(x) {
  ri <- predict_rinet(x, log_scale = TRUE)[[1]]$reference_interval
  c(unname(ri["lower"]), unname(ri["upper"]))
}))

# (2) Python re-implementation used in the study (point-estimate part of rinet_runner.py),
#     using the model and scaler files shipped inside the installed R package
py$model_dir <- system.file("models", package = "rinet")
py$x_list <- lapply(datasets, as.numeric)
py_run_string("
import numpy as np, pickle, keras
from scipy.stats import norm
model = keras.models.load_model(model_dir + '/rinet_1d.keras')
scaler = pickle.load(open(model_dir + '/scaler_1d.pkl', 'rb'))
BREAKS = np.linspace(-4, 4, 101); Z = norm.ppf([0.025, 0.975])
def features(logx):
    m = logx.mean(); s = logx.std()                      # population SD, as in the R package
    h, _ = np.histogram((logx - m) / s, bins=BREAKS, density=True)
    h = (h - h.min()) / (h.max() - h.min())
    return h, m, s
py_res = []
for x in x_list:
    lx = np.log(np.asarray(x, dtype=float))
    h, m, s = features(lx)
    pred = scaler.inverse_transform(model.predict(h[None, :, None], verbose=0))[0]
    mu = pred[0] * s + m; sd = pred[1] * s
    py_res.append([float(np.exp(mu + Z[0] * sd)), float(np.exp(mu + Z[1] * sd))])
")
p_res <- do.call(rbind, py$py_res)

# (3) Comparison
rel_LL <- abs(r_res[, 1] / p_res[, 1] - 1)
rel_UL <- abs(r_res[, 2] / p_res[, 2] - 1)
res <- data.frame(dataset = seq_along(datasets), n = sapply(datasets, length),
                  R_LL = r_res[, 1], Py_LL = p_res[, 1], R_UL = r_res[, 2], Py_UL = p_res[, 2],
                  rel_diff_LL = rel_LL, rel_diff_UL = rel_UL)
print(head(res, 10), digits = 7)
cat("\n==== RESULT ====\n")
cat("rinet version:", as.character(packageVersion("rinet")), "\n")
cat("Datasets compared:", nrow(res), "\n")
cat("Max relative difference, lower limit:", format(max(rel_LL), digits = 3), "\n")
cat("Max relative difference, upper limit:", format(max(rel_UL), digits = 3), "\n")
cat(if (max(rel_LL, rel_UL) < 1e-5) "AGREEMENT: Python re-implementation reproduces the R package.\n"
    else "DIFFERENCE FOUND: please send the full output.\n")
write.csv(res, "rinet_R_vs_python.csv", row.names = FALSE)
