# =============================================================================
# 01_nys_roc_pr_analysis.R
#
# NYS breast cancer case study (thesis, Chapter 6.1 / Appendix A).
#
# Treats the per-location relative risk (REL_RISK) estimated by SaTScan as a
# score and evaluates how well it discriminates the locations inside each of
# the five detected clusters from the rest, using:
#   * ROC curve + AUC (with automatic score direction for low-rate clusters)
#   * Youden-optimal threshold and confusion-matrix metrics at that threshold
#   * Precision-Recall curve + AUPRC
#
# Inputs  (data/nys/):
#   NYS_BreastCancer_results.rr.dbf.xlsx  SaTScan per-location RR (LOC_ID, REL_RISK)
#   NYS_BreastCancer.geo                  LOC_ID LAT LON (space separated)
#   cluster1_locations.kml ... cluster5_locations.kml  SaTScan KML output
#
# Outputs (outputs/nys/):
#   ROC_<cluster>.png, PR_<cluster>.png
#   roc_summary_all_clusters.csv
#   roc_summary_all_clusters.tex
#
# Run from the repository root:  Rscript R/01_nys_roc_pr_analysis.R
# =============================================================================


# ---- Packages ---------------------------------------------------------------
req <- c(
  "readxl", "dplyr", "stringr", "ggplot2",
  "pROC", "PRROC", "purrr", "tibble", "readr", "RANN"
)
to_install <- setdiff(req, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(pROC)
library(PRROC)
library(purrr)
library(tibble)
library(readr)
library(RANN)


# ---- Paths ------------------------------------------------------------------
data_dir   <- file.path("data", "nys")
excel_path <- file.path(data_dir, "NYS_BreastCancer_results.rr.dbf.xlsx")
geo_path   <- file.path(data_dir, "NYS_BreastCancer.geo")
kml_names  <- sprintf("cluster%d_locations.kml", 1:5)
outdir     <- file.path("outputs", "nys")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# ---- Load per-location SaTScan rr.dbf (Excel) -------------------------------
df <- read_xlsx(excel_path)
names(df) <- toupper(trimws(names(df)))

if (!all(c("LOC_ID", "REL_RISK") %in% names(df))) {
  stop("Excel must contain columns LOC_ID and REL_RISK.")
}

df <- df %>%
  mutate(
    LOC_ID   = as.character(LOC_ID) %>%
      str_replace_all("\\D", "") %>%
      na_if("") %>%
      as.numeric(),
    REL_RISK = suppressWarnings(as.numeric(REL_RISK))
  ) %>%
  filter(!is.na(LOC_ID), !is.na(REL_RISK))

known_ids <- unique(df$LOC_ID)


# ---- Load GEO file (LOC_ID, LAT, LON) ---------------------------------------
geo <- read.table(geo_path, header = FALSE, stringsAsFactors = FALSE)
names(geo) <- c("LOC_ID", "LAT", "LON")
geo$LOC_ID <- as.numeric(geo$LOC_ID)

# Coordinate matrix for nearest-neighbour search
geo_mat <- as.matrix(geo[, c("LAT", "LON")])


# ---- Extract LOC_IDs from a SaTScan KML via coordinates + nearest neighbour --
extract_ids_from_kml <- function(path, known_ids, geo, geo_mat) {
  if (!file.exists(path)) {
    message("[WARN] Missing KML: ", path)
    return(numeric(0))
  }

  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")

  # Match <coordinates>lon,lat,0</coordinates>
  coord_blocks <- str_extract_all(
    txt,
    "<coordinates>[-0-9\\.]+,[-0-9\\.]+,0</coordinates>"
  )[[1]]

  if (length(coord_blocks) == 0) {
    message("[WARN] No <coordinates> found in ", basename(path))
    return(numeric(0))
  }

  mat <- str_match(
    coord_blocks,
    "<coordinates>([-0-9\\.]+),([-0-9\\.]+),0</coordinates>"
  )
  lon <- as.numeric(mat[, 2])
  lat <- as.numeric(mat[, 3])

  # Query points in (lat, lon) order to match the GEO matrix
  query_mat <- cbind(lat, lon)

  # For each KML point find the closest GEO point
  nn      <- RANN::nn2(data = geo_mat, query = query_mat, k = 1)
  loc_ids <- geo$LOC_ID[nn$nn.idx[, 1]]

  # Keep only IDs that exist in the SaTScan rr table
  loc_ids <- unique(loc_ids[loc_ids %in% known_ids])

  message("Parsed ", length(loc_ids), " LOC_IDs from ", basename(path))
  loc_ids
}

# cluster_key -> vector of LOC_IDs
cluster_sets <- set_names(
  vector("list", length(kml_names)),
  tools::file_path_sans_ext(kml_names)
)
for (i in seq_along(kml_names)) {
  key <- tools::file_path_sans_ext(kml_names[i])  # e.g. "cluster1_locations"
  cluster_sets[[key]] <- extract_ids_from_kml(
    file.path(data_dir, kml_names[i]), known_ids, geo, geo_mat
  )
}


# ---- Helpers: ROC / PR and Youden index -------------------------------------

# Fit ROC on REL_RISK and on -REL_RISK (low-rate clusters), keep the better one
compute_roc_dir <- function(scores, labels) {
  roc1 <- roc(labels,  scores, quiet = TRUE, direction = ">")
  roc2 <- roc(labels, -scores, quiet = TRUE, direction = ">")
  if (as.numeric(auc(roc2)) > as.numeric(auc(roc1))) {
    list(roc = roc2, direction = "-REL_RISK", use_scores = -scores)
  } else {
    list(roc = roc1, direction = "REL_RISK",  use_scores =  scores)
  }
}

youden_stats <- function(roc_obj) {
  best <- coords(
    roc_obj, "best", best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  thr  <- best$threshold
  y    <- roc_obj$response
  s    <- roc_obj$predictor
  pred <- ifelse(s >= thr, 1, 0)
  tp <- sum(pred == 1 & y == 1)
  fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1)
  tn <- sum(pred == 0 & y == 0)
  list(
    tau  = thr,
    TPR  = as.numeric(best$sensitivity),
    FPR  = as.numeric(1 - best$specificity),
    SPEC = as.numeric(best$specificity),
    PREC = ifelse(tp + fp > 0, tp / (tp + fp), NA_real_),
    TP = tp, FP = fp, FN = fn, TN = tn
  )
}

compute_pr <- function(scores, labels) {
  sc_pos <- scores[labels == 1]
  sc_neg <- scores[labels == 0]
  if (length(sc_pos) == 0 || length(sc_neg) == 0) {
    return(list(curve = NULL, auprc = NA_real_))
  }
  pr <- PRROC::pr.curve(
    scores.class0 = sc_pos,
    scores.class1 = sc_neg,
    curve = TRUE
  )
  list(curve = pr$curve, auprc = as.numeric(pr$auc.integral))
}


# ---- Evaluate each cluster --------------------------------------------------
sum_rows <- list()

for (key in names(cluster_sets)) {
  pos_ids <- cluster_sets[[key]]
  labels  <- ifelse(df$LOC_ID %in% pos_ids, 1, 0)
  scores  <- df$REL_RISK

  # Need at least some positives and negatives
  if (length(unique(labels)) < 2) {
    message("[WARN] Cluster ", key,
            " has <2 label levels (maybe 0 LOC_IDs matched). Skipping.")
    next
  }

  rocdir     <- compute_roc_dir(scores, labels)
  roc_obj    <- rocdir$roc
  direction  <- rocdir$direction
  use_scores <- rocdir$use_scores

  auc_val <- as.numeric(auc(roc_obj))
  ystats  <- youden_stats(roc_obj)
  pr      <- compute_pr(use_scores, labels)

  # --- ROC plot ---
  roc_df <- tibble(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  )
  p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey60") +
    geom_line(linewidth = 1) +
    labs(
      title    = paste0("ROC - ", key, " (", direction, ")"),
      subtitle = sprintf("AUC = %.3f", auc_val),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal(base_size = 11)
  roc_path <- file.path(outdir, paste0("ROC_", key, ".png"))
  ggsave(roc_path, p_roc, width = 5, height = 5, dpi = 180)

  # --- PR plot ---
  if (!is.null(pr$curve)) {
    pr_df <- tibble(
      Recall    = pr$curve[, 1],
      Precision = pr$curve[, 2]
    )
    p_pr <- ggplot(pr_df, aes(x = Recall, y = Precision)) +
      geom_line(linewidth = 1) +
      labs(
        title    = paste0("Precision-Recall - ", key, " (", direction, ")"),
        subtitle = sprintf("AUPRC = %.3f", pr$auprc),
        x = "Recall (TPR)",
        y = "Precision"
      ) +
      theme_minimal(base_size = 11)
  } else {
    p_pr <- ggplot() + theme_void() +
      ggtitle(paste0("PR curve not available - ", key))
  }
  pr_path <- file.path(outdir, paste0("PR_", key, ".png"))
  ggsave(pr_path, p_pr, width = 5, height = 5, dpi = 180)

  # --- Collect metrics ---
  sum_rows[[length(sum_rows) + 1]] <- tibble(
    cluster_key        = key,
    n_total            = length(labels),
    n_pos              = sum(labels == 1),
    n_neg              = sum(labels == 0),
    score_direction    = direction,
    AUC                = round(auc_val, 4),
    AUPRC              = round(pr$auprc, 4),
    tau_Youden         = ystats$tau,
    TPR_at_tau         = round(ystats$TPR, 4),
    FPR_at_tau         = round(ystats$FPR, 4),
    Specificity_at_tau = round(ystats$SPEC, 4),
    Precision_at_tau   = round(ystats$PREC, 4),
    TP                 = ystats$TP,
    FP                 = ystats$FP,
    FN                 = ystats$FN,
    TN                 = ystats$TN,
    roc_png            = roc_path,
    pr_png             = pr_path,
    n_ids_in_kml       = length(pos_ids)
  )
}

summary_tbl <- bind_rows(sum_rows)


# ---- Save CSV + LaTeX summary -----------------------------------------------
csv_path <- file.path(outdir, "roc_summary_all_clusters.csv")
write_csv(summary_tbl, csv_path)

if (nrow(summary_tbl) > 0) {
  tex_lines <- c(
    "\\begin{tabular}{l r r r l r r r r r}",
    "\\toprule",
    paste0("Cluster & $N$ & Pos & Neg & Score & AUC & AUPRC & ",
           "TPR$_{\\ast}$ & FPR$_{\\ast}$ & $\\tau_{\\ast}$\\\\"),
    "\\midrule"
  )
  for (i in seq_len(nrow(summary_tbl))) {
    r <- summary_tbl[i, ]
    tex_lines <- c(
      tex_lines,
      sprintf(
        "%s & %d & %d & %d & %s & %.3f & %.3f & %.3f & %.3f & %s\\\\",
        r$cluster_key, r$n_total, r$n_pos, r$n_neg, r$score_direction,
        r$AUC, r$AUPRC, r$TPR_at_tau, r$FPR_at_tau,
        ifelse(is.na(r$tau_Youden), "-", format(r$tau_Youden, digits = 5))
      )
    )
  }
  tex_lines <- c(tex_lines, "\\bottomrule", "\\end{tabular}")
  tex_path  <- file.path(outdir, "roc_summary_all_clusters.tex")
  writeLines(tex_lines, tex_path)

  message("Done.")
  message("Summary CSV: ", csv_path)
  message("LaTeX table: ", tex_path)
  for (i in seq_len(nrow(summary_tbl))) {
    cat(sprintf("- %s: ROC=%s | PR=%s\n",
                summary_tbl$cluster_key[i],
                summary_tbl$roc_png[i],
                summary_tbl$pr_png[i]))
  }
} else {
  message("No clusters with valid labels found; summary table not created.")
}
