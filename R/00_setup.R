# Paths to the raw PhysioNet files (csv.gz) and output folders.
# MIMIC_DIR: mimic-iv-3.1 folder (contains hosp/ and icu/)
MIMIC_DIR <- Sys.getenv("MIMIC_DIR", "data/mimic-iv-3.1")
WORK_DIR  <- "work"
OUT_DIR   <- "results"
dir.create(WORK_DIR, showWarnings = FALSE)
dir.create(OUT_DIR, showWarnings = FALSE)
pkgs <- c("DBI", "duckdb", "data.table", "survival", "cmprsk", "splines", "mice", "pROC",
          "EValue", "WeightIt", "MatchIt", "ggplot2", "patchwork", "scales", "jsonlite")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install packages: ", paste(miss, collapse = ", "))
suppressMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
