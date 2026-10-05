# =============================================================================
# 02_simulation_spatial_overlap.R
#
# Simulation study on a 20 x 20 lattice (thesis, Chapter 6.2 / Appendix B).
#
# Compares the regions detected by SaTScan with the known (true) high-rate and
# low-rate clusters using spatial overlap metrics:
#   sensitivity, precision, Jaccard index and symmetric difference area.
#
# Inputs  (data/simulation/):
#   sim_example_true_labels.csv     ";"-separated: loc_id, in_high_true, in_low_true
#   sim_example_results.gis.dbf     SaTScan GIS output (LOC_ID, CLUSTER, LOC_RR)
#
# Output  (outputs/simulation/):
#   sim_overlap_summary.csv
#
# Run from the repository root:  Rscript R/02_simulation_spatial_overlap.R
# =============================================================================


# ---- Packages ---------------------------------------------------------------
req <- c("dplyr", "readr", "foreign", "tibble")
to_install <- setdiff(req, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

library(dplyr)
library(readr)
library(foreign)   # read.dbf
library(tibble)


# ---- Paths ------------------------------------------------------------------
data_dir <- file.path("data", "simulation")
outdir   <- file.path("outputs", "simulation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# ---- 1. Ground-truth labels -------------------------------------------------
true_labels <- read_delim(
  file.path(data_dir, "sim_example_true_labels.csv"),
  delim = ";",
  show_col_types = FALSE,
  trim_ws = TRUE
) %>%
  mutate(
    loc_id          = as.character(loc_id),
    in_high_true    = as.logical(as.integer(in_high_true)),
    in_low_true     = as.logical(as.integer(in_low_true)),
    in_cluster_true = in_high_true | in_low_true   # any true cluster
  )


# ---- 2. SaTScan GIS output --------------------------------------------------
gis_out <- read.dbf(file.path(data_dir, "sim_example_results.gis.dbf"),
                    as.is = TRUE) %>%
  as_tibble() %>%
  mutate(
    LOC_ID  = as.character(LOC_ID),
    CLUSTER = as.integer(CLUSTER),
    LOC_RR  = as.numeric(LOC_RR)
  )


# ---- 3. Merge truth + detections by location ID -----------------------------
sim_merged <- true_labels %>%
  inner_join(
    gis_out %>% select(LOC_ID, CLUSTER, LOC_RR),
    by = c("loc_id" = "LOC_ID")
  ) %>%
  mutate(
    in_scan_any  = CLUSTER > 0,                  # any detected cluster
    in_scan_high = (CLUSTER > 0) & (LOC_RR > 1), # detected high-rate region
    in_scan_low  = (CLUSTER > 0) & (LOC_RR < 1)  # detected low-rate region
  )


# ---- 4. Overlap metrics -----------------------------------------------------
overlap_summary <- function(df, true_col, pred_col, label) {
  y_true <- df[[true_col]]
  y_pred <- df[[pred_col]]

  TP <- sum( y_true &  y_pred)
  FP <- sum(!y_true &  y_pred)
  FN <- sum( y_true & !y_pred)
  TN <- sum(!y_true & !y_pred)

  tibble(
    cluster_type = label,
    TP = TP, FP = FP, FN = FN, TN = TN,
    sensitivity = TP / (TP + FN),
    precision   = TP / (TP + FP),
    jaccard     = TP / (TP + FP + FN),
    symdiff     = FP + FN
  )
}

sim_overlap_summary <- bind_rows(
  overlap_summary(sim_merged, "in_high_true",    "in_scan_high", "High-rate"),
  overlap_summary(sim_merged, "in_low_true",     "in_scan_low",  "Low-rate"),
  overlap_summary(sim_merged, "in_cluster_true", "in_scan_any",  "Any cluster")
)


# ---- 5. Save and print ------------------------------------------------------
write_csv(sim_overlap_summary, file.path(outdir, "sim_overlap_summary.csv"))
print(sim_overlap_summary)
