# Occupancy-Based Building Clustering =========================================
# Description:
#   This script clusters buildings by their occupancy patterns
#   using three methods:
#     1. CFDA  – Categorical Functional Data Analysis (Preda et al., 2021)
#     2. ClickClust – Model-based clustering of categorical sequences
#                     (Melnykov, 2016)
#     3. FunLBM – Functional Latent Block Model for co-clustering of buildings
#                 and time periods (Bouveyron et al., 2018)
#
#   Cluster solutions are compared via cross-tabulation and the Adjusted Rand
#   Index (ARI). All figures are saved to plots/ as PDF files.
#
# Input:
#   data/01_large_scale_occ.fst  –  Hourly binary occupancy estimates for all
#                                   buildings (columns: bldg_id, time,
#                                   occ_estimated)
# Output:
#   plots/01_hourly_occupancy.pdf        –  Annual occupancy heatmap (Fig. 7)
#   plots/02_hourly_occupancy_cfda.pdf   –  CFDA cluster heatmaps
#   plots/03_hourly_occupancy_lbm.pdf    –  FunLBM cluster heatmaps
#   plots/04_lbm_combined.pdf            –  FunLBM time clusters + profiles
#   plots/05_cls_comparison.pdf          –  Cross-method cluster comparison
#
# References:
#   Preda et al. (2021) https://doi.org/10.3390/math9233074
#   Melnykov (2016)     https://doi.org/10.18637/jss.v074.i09
#   Bouveyron et al. (2018) https://doi.org/10.1111/rssc.12260
#
# Author:  Markus Schaffer (msch@build.aau.dk)


# 1. Setup =====================================================================


## 1.1 Packages ---------------------------------------------------------------

library(ClickClust)
library(cfda)
library(funLBM)
library(fda)
library(NbClust)
library(plyr)
library(data.table)
library(fst)
library(purrr)
library(lubridate)
library(furrr)
library(progressr)
library(ggplot2)
library(patchwork)
library(ggdendro)
library(scales)
library(ggpattern)
library(mclust, include.only = "adjustedRandIndex")


## 1.2 Output directories -----------------------------------------------------

if (!dir.exists("data/cluster")) {
  dir.create("data/cluster", recursive = TRUE)
}
if (!dir.exists("data/cluster/tmp")) {
  dir.create("data/cluster/tmp", recursive = TRUE)
}
if (!dir.exists("plots")) {
  dir.create("plots", recursive = TRUE)
}


## 1.3 Plot list (all figures accumulated here) --------------------------------

plots <- list()


# 2. Helper definitions ========================================================

## 2.1 School / public holidays (Denmark, 2022) --------------------------------
#
# Returns a character label for each POSIXct timestamp.  Dates that fall
# outside any holiday window are labelled "No holiday".

holidays <- function(date) {
  fcase(
    date %within% interval(
      as.POSIXct("2022-01-01", tz = "Europe/Copenhagen"),
      as.POSIXct("2022-01-03 23:59", tz = "Europe/Copenhagen")
    ) |
      date %within% interval(
        as.POSIXct("2022-12-22", tz = "Europe/Copenhagen"),
        as.POSIXct("2023-01-01 23:59", tz = "Europe/Copenhagen")
      ),
    "Christmas holiday",
    date %within% interval(
      as.POSIXct("2022-02-19", tz = "Europe/Copenhagen"),
      as.POSIXct("2022-02-27 23:59", tz = "Europe/Copenhagen")
    ),
    "Winter holiday",
    date %within% interval(
      as.POSIXct("2022-04-09", tz = "Europe/Copenhagen"),
      as.POSIXct("2022-04-18 23:59", tz = "Europe/Copenhagen")
    ),
    "Easter holiday",
    date %within% interval(
      as.POSIXct("2022-06-25", tz = "Europe/Copenhagen"),
      as.POSIXct("2022-08-07 23:59", tz = "Europe/Copenhagen")
    ),
    "Summer holiday",
    date %within% interval(
      as.POSIXct("2022-10-15", tz = "Europe/Copenhagen"),
      as.POSIXct("2022-10-23 23:59", tz = "Europe/Copenhagen")
    ),
    "Autumn holiday",
    date %within% interval(
      as.POSIXct("2022-05-26", tz = "Europe/Copenhagen"),
      as.POSIXct("2022-05-29 23:59", tz = "Europe/Copenhagen")
    ),
    "\nAscension holiday",
    date %within% interval(
      as.POSIXct("2022-05-13", tz = "Europe/Copenhagen"),
      as.POSIXct("2022-05-15 23:59", tz = "Europe/Copenhagen")
    ),
    "Great prayer holiday",
    date %within% interval(
      as.POSIXct("2022-06-04", tz = "Europe/Copenhagen"),
      as.POSIXct("2022-06-06 23:59", tz = "Europe/Copenhagen")
    ),
    "\n\nWhitsunday holiday",
    default = "No holiday"
  )
}


## 2.2 Colour palettes --------------------------------------------------------

gradient_occ <- c(
  "#F5F3C1", "#EAF0B5", "#DDECBF", "#D0E7CA", "#C2E3D2", "#B5DDD8",
  "#A8D8DC", "#9BD2E1", "#8DCBE4", "#81C4E7", "#7BBCE7", "#7EB2E4",
  "#88A5DD", "#9398D2", "#9B8AC4", "#9D7DB2", "#9A709E", "#906388",
  "#805770", "#684957", "#46353A"
)

# Discrete colours for FunLBM time clusters (T1–T5)
lbm_temp_colors <- c("#0077BB", "#33BBEE", "#EE7733", "#009988", "#CC3311")


## 2.3 ggplot2 themes ---------------------------------------------------------

# Shared formatting applied by both custom themes below
.theme_common <- function() {
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    line = element_line(colour = "black", linewidth = 0.1),
    rect = element_rect(colour = "black", linewidth = 0.1),
    strip.background = element_rect(fill = "white", linewidth = 0.1),
    legend.position = "bottom"
  )
}

# Clean black/white theme – used for static comparison and selection plots
theme_nice <- function() {
  theme_bw() +
    .theme_common() +
    theme(
      axis.ticks   = element_line(colour = "black", linewidth = 0.1),
      panel.border = element_rect(colour = "black", linewidth = 0.1)
    )
}

# Classic theme – used for time-series and faceted heatmap plots
theme_time <- function() {
  theme_classic() +
    .theme_common() +
    theme(
      panel.spacing.x = unit(0, "mm"),
      panel.spacing.y = unit(1.5, "mm"),
      legend.title = element_text(size = 6, colour = "black"),
      axis.ticks = element_line(colour = "black", linewidth = 0.1),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(size = 4),
      legend.position = "bottom"
    )
}

# 3. Load and preprocess data =================================================

occ_dt <- read_fst("data/01_large_scale_occ.fst", as.data.table = TRUE)

# Remove the first and last calendar days, which are incomplete in the dataset
occ_dt[, day := yday(time)]
occ_dt <- occ_dt[day != 1 & day != 365]


## 3.1 Daylight-saving-time correction ----------------------------------------
#
# During the DST transition, clocks skip one hour, producing a 23-hour day.
# To keep all days at exactly 24 observations per building, the hour just
# before the DST end is duplicated and the DST end hour is removed.  A
# monotone "pseudo_hour" index (1–24) is then assigned so that downstream
# analyses always see 24 evenly-spaced time steps per day.

get_dst <- function(y = 2019, tz) {
  start <- map_chr(y, ~ paste0(.x, "-01-01"))
  end <- map_chr(y, ~ paste0(.x, "-12-31"))
  d1 <- map2(start, end, ~ seq(as.POSIXct(.x, tz = tz),
    as.POSIXct(.y, tz = tz),
    by = "hour"
  ))
  map_dfr(d1, ~ range(.x[dst(.x)]) |> setNames(c("start", "end")))
}

dst_range <- get_dst(2022, "Europe/Copenhagen")

occ_dt <- rbind(occ_dt, occ_dt[time == dst_range$start]) # duplicate hour before DST end
occ_dt <- occ_dt[time != dst_range$end] # drop the DST end hour
setorder(occ_dt, bldg_id, time)

# Assign sequential pseudo-hours; override the DST day to ensure 1–24
occ_dt[, pseudo_hour := hour(time)]
occ_dt[day == yday(dst_range$start), pseudo_hour := seq_len(.N) - 1, by = bldg_id]


## 3.2 Calendar variables -----------------------------------------------------

occ_dt[, month := month(time,
  label = TRUE, abbr = TRUE,
  locale = "English_United States"
)]
occ_dt[, weekday := wday(time,
  label = TRUE, week_start = 1,
  locale = "English_United States"
)]


# 4. Exploratory visualisation – annual occupancy overview ====================

## 4.1 Holiday annotation strip -----------------------------------------------

holiday_dt <- occ_dt[bldg_id == bldg_id[1], .(time, month, day)]
holiday_dt[, holiday := holidays(time)]

# Compute x-extent of each holiday period; split Christmas across Jan/Dec
holiday_dt <- rbind(
  holiday_dt[holiday != "No holiday" & holiday != "Christmas holiday",
    .(x_min = min(day), x_max = max(day)),
    by = holiday
  ],
  holiday_dt[holiday == "Christmas holiday" & month == "Jan",
    .(x_min = min(day), x_max = max(day)),
    by = holiday
  ],
  holiday_dt[holiday == "Christmas holiday" & month == "Dec",
    .(x_min = min(day), x_max = max(day)),
    by = holiday
  ]
)

holiday_dt[, x_min := x_min - 0.5]
holiday_dt[, x_max := x_max + 0.5]
holiday_dt[, y_min := 0]
holiday_dt[, y_max := 1]
holiday_dt[, breaks := x_min + (x_max - x_min) / 2]
holiday_dt[, label := gsub(" holiday", "", holiday)]

plots[["holidays"]] <- ggplot() +
  geom_rect(
    data = holiday_dt,
    aes(
      xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max,
      fill = I("transparent"), colour = I("black")
    ),
    linewidth = 0.15
  ) +
  scale_y_discrete(name = NULL, expand = c(0, 0)) +
  scale_x_continuous(
    breaks       = holiday_dt$breaks,
    minor_breaks = NULL,
    labels       = holiday_dt$label,
    expand       = c(0, 0),
    name         = NULL,
    guide        = guide_axis(angle = 0)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 5, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.1, linetype = 1),
    panel.grid.major = element_blank()
  )


## 4.2 Annual hourly-occupancy heatmap (all buildings combined) ----------------

hourly_occ <- occ_dt[, .(mean_occ = mean(occ_estimated)), by = .(month, day, pseudo_hour)]

plots[["hourly_occ"]] <- ggplot() +
  geom_raster(
    data = hourly_occ,
    aes(x = day, y = pseudo_hour, fill = mean_occ), vjust = 0
  ) +
  scale_x_continuous(breaks = NULL, expand = c(0, 0), name = NULL) +
  scale_y_reverse(name = "hour", expand = c(0, 0), breaks = seq(0, 24)) +
  facet_grid(cols = vars(month), scale = "free", space = "free") +
  scale_fill_gradientn(
    colours = rev(gradient_occ),
    name = "Buildings occupied",
    na.value = "transparent",
    limits = c(0.3, 1),
    oob = scales::squish,
    breaks = c(0.3, 0.5, 0.7, 0.9),
    labels = paste0(c("\U2264 ", rep("", 3)), c(30, 50, 70, 90), "%")
  ) +
  scale_color_manual(values = "black", label = "holidays", name = NULL) +
  guides(
    color = guide_legend(override.aes = list(fill = "white", key.size = unit(4, "mm"))),
    fill  = guide_colorbar(title.position = "bottom", title.hjust = 0.5)
  ) +
  theme_time()


## 4.3 Overlay heatmap and holiday strip (shared layout) ----------------------

shared_layout <- c(area(1, 1), area(1, 1))

plots[["hourly_occ_overview"]] <- plots[["hourly_occ"]] +
  plots[["holidays"]] +
  plot_layout(design = shared_layout)

# Note: the pdf device produces better fill alignment than cairo_pdf but cannot
# render the ≤ Unicode character.  The published figure used pdf() with the ≤
# inserted manually in post-processing.
# This figure is not shown in the manuscript due to space issues
ggsave(
  "plots/01_hourly_occupancy.pdf",
  plot = plots[["hourly_occ_overview"]],
  width = 122.5,
  height = 80,
  units = "mm",
  device = pdf
)


# 5. Clustering method 1: CFDA ================================================
#
# Reference:
#   Preda, C., Grimonprez, Q., & Vandewalle, V. (2021).
#   Categorical functional data analysis. The cfda R package.
#   Mathematics, 9(23), 1–31. https://doi.org/10.3390/math9233074


## 5.1 Reformat data for cfda -------------------------------------------------
#
# cfda expects a long-format data.table with columns (id, time, state), where
# `time` is a consecutive integer index within each building.

cfda_dt <- occ_dt[, .(id = bldg_id, state = occ_estimated)]
cfda_dt[, time := seq_len(.N), by = id]


## 5.2 Compute optimal encoding -----------------------------------------------
#
# A B-spline basis with nbasis = 10 is used to represent the occupancy
# trajectories.  Sensitivity analyses in the original paper showed that varying
# nbasis between 5, 10, and 20 has negligible effect on the resulting encoding.

cfda_basis <- create.bspline.basis(
  rangeval = c(1, cfda_dt[id == id[1], max(time)]),
  nbasis   = 10,
  norder   = 4
)

set.seed(123)
cfda_encoding <- compute_optimal_encoding(
  cfda_dt, cfda_basis,
  computeCI = FALSE,
  method = "parallel",
  nCores = 16
)
saveRDS(cfda_encoding, "data/cluster/cfda.RDS")
cfda_encoding <- readRDS("data/cluster/cfda.RDS")


## 5.3 Select number of principal components and clusters ---------------------

# Retain the minimum number of PCs that explain ≥ 90 % of variance
plot(cumsum(prop.table(cfda_encoding$eigenvalues)))
cfda_nPc90 <- which(cumsum(prop.table(cfda_encoding$eigenvalues)) > 0.90)[1]

# Determine optimal cluster count via the Calinski–Harabasz (CH) index,
# testing k = 2 … 15 with Ward's D² linkage
set.seed(123)
cfda_result <- NbClust(
  cfda_encoding$pc[, 1:cfda_nPc90],
  distance = "euclidean",
  method   = "ward.D2",
  max.nc   = 15,
  index    = "ch"
)
cfda_result$Best.nc # CH index selects k = 2


# Manual selection via dendrogram

cfda_tree <- hclust(dist(cfda_encoding$pc[, 1:cfda_nPc90]), method = "ward.D2")

plots$cfda_dendo <- ggplot(
  segment(dendro_data(as.dendrogram(cfda_tree), type = "rectangle"))
) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_nice()

plots$cfda_dendo
# → clusters = 4


cfda_res <- cutree(cfda_tree, k = 4) |>
  as.data.table(keep.rownames = TRUE)
setnames(cfda_res, c("V1", "V2"), c("bldg_id", "cfda"))

# Clusters are relabelled O1–O4 in descending order of size so that O1 always
# refers to the largest group.
cfda_res[, bldg_id := as.numeric(bldg_id)]
cfda_res[, cfda_n := .N, by = cfda]
setorder(cfda_res, -cfda_n)
cfda_res[, cfda := paste0("O", rleid(cfda_n))]

occ_dt <- merge(occ_dt, cfda_res, by = "bldg_id")


## 5.6 CFDA cluster heatmap ---------------------------------------------------

hourly_occ_cfda <- occ_dt[, .(mean_occ = mean(occ_estimated)),
  by = .(month, day, pseudo_hour, cfda, cfda_n)
]
hourly_occ_cfda[, cfda := paste0(cfda, "\nn = ", cfda_n)]

plots[["hourly_occ_cfda"]] <- ggplot() +
  geom_raster(
    data = hourly_occ_cfda,
    aes(x = day, y = pseudo_hour, fill = mean_occ), vjust = 0
  ) +
  scale_x_continuous(breaks = NULL, expand = c(0, 0), name = NULL) +
  scale_y_reverse(name = "hour", expand = c(0, 0), breaks = c(0, 12, 24)) +
  facet_grid(cols = vars(month), rows = vars(cfda), scale = "free", space = "free") +
  scale_fill_gradientn(
    colours  = rev(gradient_occ),
    name     = "Buildings occupied",
    na.value = "transparent",
    limits   = c(0.3, 1),
    oob      = scales::squish,
    breaks   = c(0.3, 0.5, 0.7, 0.9),
    labels   = paste0(c("\U2264 ", rep("", 3)), c(30, 50, 70, 90), "%")
  ) +
  scale_color_manual(values = "black", label = "holidays", name = NULL) +
  guides(
    color = guide_legend(override.aes = list(fill = "white", key.size = unit(4, "mm"))),
    fill  = guide_colorbar(title.position = "bottom", title.hjust = 0.5)
  ) +
  theme_time()

plots[["hourly_occ_overview_cfda"]] <- plots[["hourly_occ_cfda"]] +
  plots[["holidays"]] +
  plot_layout(design = shared_layout)

ggsave(
  "plots/02_hourly_occupancy_cfda.pdf",
  plot = plots[["hourly_occ_overview_cfda"]],
  width = 122.5,
  height = 100,
  units = "mm",
  device = pdf
)


# 6. Clustering method 2: ClickClust ==========================================
#
# Reference:
#   Melnykov, V. (2016).
#   ClickClust: An R package for model-based clustering of categorical
#   sequences. Journal of Statistical Software, 74(9).
#   https://doi.org/10.18637/jss.v074.i09


## 6.1 Reformat data for ClickClust -------------------------------------------
#
# ClickClust requires a list of integer-coded state sequences (1-indexed).

clickclust_lst <- split(occ_dt[, .(occ_estimated, bldg_id)],
  by = "bldg_id", keep.by = FALSE
)
clickclust_lst <- map(clickclust_lst, ~ as.numeric(.x[, occ_estimated]) + 1L)
clickclust_obj <- click.read(clickclust_lst)


## 6.2 Fit ClickClust models for k = 2 … 20 -----------------------------------
#
# BIC is used for model selection.

plan(multisession, workers = 8)
clickclust_cluster <- future_map(
  2:20,
  ~ click.EM(
    X = clickclust_obj$X, K = .x, iter = 10, r = 1000,
    scale.const = 1
  ),
  .options = furrr_options(chunk_size = 1, seed = 123)
)
plan(sequential)

saveRDS(clickclust_cluster, "data/cluster/clickclust.RDS")
clickclust_cluster <- readRDS("data/cluster/clickclust.RDS")


## 6.3 Summarise model selection criteria ------------------------------------

clickclust_res <- map(
  clickclust_cluster, ~ data.table(
    loglik = .x$logl,
    bic = .x$BIC,
    k = length(.x$alpha)
  )
) |> rbindlist()

# All fail
clickclust_res


# 7. Clustering method 3: FunLBM ==============================================
#
# Reference:
#   Bouveyron, C., Bozzi, L., Jacques, J., & Jollois, F.-X. (2018).
#   The functional latent block model for the co-clustering of electricity
#   consumption curves. Journal of the Royal Statistical Society, Series C,
#   67(4), 897–915. https://doi.org/10.1111/rssc.12260


## 7.1 Reformat data for FunLBM -----------------------------------------------
#
# FunLBM expects a 3-D array of dimensions [buildings × days × hours].

lbm_dt <- occ_dt[, .(occ_estimated = as.numeric(occ_estimated), bldg_id)]

timepoints <- 24L
bldgs <- lbm_dt[, uniqueN(bldg_id)]
days <- lbm_dt[, .N] / timepoints / bldgs

lbm_occ_array <- array(lbm_dt$occ_estimated, dim = c(timepoints, days, bldgs))
lbm_occ_array <- aperm(lbm_occ_array, c(3, 2, 1)) # → [buildings × days × hours]
dimensions <- dim(lbm_occ_array)


## 7.2 Select the number of Fourier basis functions ---------------------------
#
# GCV (generalised cross-validation) is computed for odd values of nbasis
# from 3 to 23 and averaged across all building–day combinations.

bspline_basis <- function(data, timepoints, nbasis) {
  basisobj <- create.fourier.basis(
    rangeval = range(seq_len(timepoints)),
    nbasis = nbasis
  )
  ys <- smooth.basis(
    argvals = seq_len(timepoints), y = data,
    fdParobj = basisobj
  )
  data.table(gcv = ys$gcv, nbasis = nbasis)
}

plan(multisession, workers = 5)
lbm_basis_result <- future_map(
  seq(3, timepoints - 1, 2), ~ bspline_basis(
    data = t(apply(lbm_occ_array, 3, cbind)),
    timepoints = timepoints,
    nbasis = .x
  ),
  .options = furrr_options(scheduling = 1L, seed = 123)
) |> rbindlist()
plan(sequential)

# Average GCV across building-days; pick nbasis at the minimum
basis_result_mean <- lbm_basis_result[, .(mean_gcv = mean(gcv)), by = nbasis]

plots[["no_basis"]] <- ggplot(basis_result_mean, aes(x = nbasis, y = mean_gcv)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  geom_point(data = basis_result_mean[mean_gcv == min(mean_gcv)], size = 2) +
  scale_x_continuous(
    breaks       = seq(3, timepoints - 1, 4),
    minor_breaks = seq(5, timepoints - 1, 4),
    name         = "no. of basis functions"
  ) +
  scale_y_continuous(name = "GCV") +
  coord_cartesian(xlim = c(3, 19)) +
  theme_nice()
plots[["no_basis"]]
# → nbasis = 7 is used in all subsequent FunLBM fits


## 7.3 Fit FunLBM over all row × column cluster combinations -----------------
#
# Row clusters (K) correspond to buildings; column clusters (L) correspond to
# time periods (days).  The wrapper catches errors for combinations that fail
# to converge and saves partial results to disk for fault-tolerance.

fn_best_cluster <- function(input_data, K = 5, L = 5, nbasis = 7) {
  lbl <- tryCatch(
    funLBM(input_data,
      K = K, L = L, nbasis = nbasis, maxit = 300, burn = 50,
      init = "funFEM", basis.name = "fourier",
      nbinit = 3, gibbs.it = 3
    ),
    error = function(e) "error"
  )
  # Checkpoint: save individual result in case of cluster failure
  saveRDS(
    data.table(lbl = list(lbl), K = K, L = L, nbasis = nbasis),
    paste0("data/cluster/tmp/LMB_", K, "_", L, "_.RDS")
  )
  pg()
  data.table(lbl = list(lbl), K = K, L = L, nbasis = nbasis)
}

lbm_param <- expand.grid(list(row_k = 2:12, col_l = 2:9))

plan(multisession, workers = 8)
with_progress({
  pg <- progressor(steps = nrow(lbm_param))
  lbm_result <- future_map2(
    lbm_param$row_k, lbm_param$col_l,
    ~ fn_best_cluster(lbm_occ_array, K = .x, L = .y, nbasis = 7),
    .options = furrr_options(seed = 123, chunk_size = 1L)
  ) |> rbindlist()
})
plan(sequential)

saveRDS(lbm_result, "data/cluster/lbm.RDS")
lbm_result <- readRDS("data/cluster/lbm.RDS")


## 7.4 Model selection via ICL ------------------------------------------------

# Extract ICL only for converged fits
lbm_result[map_lgl(lbl, ~ .x[[1]] != "error"), icl := map_dbl(lbl, ~ .x$icl)]

# Factor labels for plotting (O = building clusters, T = time clusters)
lbm_result[, L := factor(paste0("T", L), levels = paste0("T", sort(unique(L))))]
lbm_result[, K := factor(paste0("O", K), levels = paste0("O", sort(unique(K))))]

# Retrieve the best-fitting model
lbm_best_cls <- lbm_result[icl == max(icl, na.rm = TRUE), lbl][[1]]


## 7.5 Evaluate convergence ---------------------------------------------------

lbm_lst <- lbm_result[map_lgl(lbl, ~ .x[[1]] != "error"), lbl]
loglike_plot <- map(lbm_lst, ~ data.table(
  loglik = .x$loglik,
  K      = .x$K,
  L      = .x$L,
  it     = seq_along(.x$loglik)
)) |> rbindlist()
loglike_plot[, K_L := paste0(K, "_", L)]

plots[["lbm_icl_con"]] <- ggplot() +
  geom_line(
    data = loglike_plot[K != lbm_best_cls$K & L != lbm_best_cls$L],
    aes(x = it, y = loglik, group = K_L),
    linewidth = 0.2
  ) +
  geom_line(
    data = loglike_plot[K == lbm_best_cls$K & L == lbm_best_cls$L],
    aes(x = it, y = loglik, group = K_L),
    linewidth = 0.5, colour = "#DD3D2D"
  ) +
  scale_y_continuous(labels = label_scientific(digits = 3)) +
  labs(x = "iterations", y = "Complete log\u2212likelihood") +
  theme_nice()
plots[["lbm_icl_con"]]
# → converged

## 7.6 ICL selection heatmap --------------------------------------------------

plots[["lbm_icl_sel"]] <- ggplot(lbm_result, aes(x = L, y = K, fill = icl)) +
  geom_tile() +
  geom_text(
    data = lbm_result[icl == max(icl, na.rm = TRUE)],
    aes(x = L, y = K, label = "max", color = "NA"),
    size = 6 / .pt
  ) +
  scale_x_discrete(name = "time cluster", expand = c(0, 0)) +
  scale_y_discrete(name = "energy use cluster", expand = c(0, 0)) +
  scale_fill_gradientn(
    labels    = label_scientific(digits = 3),
    colors    = rev(gradient_occ),
    na.value  = "#FFFFFF",
    name      = "ICL"
  ) +
  scale_color_manual(values = "black", name = "Failed", labels = NULL) +
  guides(color = guide_legend(
    override.aes = list(fill = "white", color = "white"), order = 2
  )) +
  theme_nice() +
  theme(
    legend.position      = "right",
    legend.key           = element_rect(colour = "black"),
    legend.box.margin    = margin(0, 0, -15, 0)
  )
plots[["lbm_icl_sel"]]
# → 6 x 5 clusters

## 7.7 Assign cluster labels to observations ----------------------------------
#
# Column (time-period) clusters are manually relabelled T1–T5 to match the
# temporal ordering described in the paper.  Row (building) clusters are
# relabelled O1–O6 in descending order of size.

lbm_col_clust <- data.table(col = lbm_best_cls$col_clust)
lbm_col_clust[, new_col := fcase(
  col == 1, "T2",
  col == 2, "T1",
  col == 3, "T3",
  col == 4, "T4",
  col == 5, "T5"
)]

lbm_row_clust <- data.table(row = lbm_best_cls$row_clust)
lbm_row_clust[, idx := .I]
lbm_row_clust[, N := .N, by = row]
lbm_row_clust_n <- unique(lbm_row_clust[, .(row, N)])
setorder(lbm_row_clust_n, -N)
lbm_row_clust_n[, new_row := paste0("O", frank(-N))]
lbm_row_clust <- merge(lbm_row_clust, lbm_row_clust_n[, .(row, new_row)],
  by = "row", sort = FALSE
)

# Attach cluster assignments to the main data.table
occ_dt[, lbm_period := rep(lbm_col_clust$new_col, each = 24), by = bldg_id]
occ_dt[, lbm_cluster := rep(
  lbm_row_clust$new_row,
  occ_dt[, .N, by = bldg_id][, N]
)]


## 7.8 FunLBM cluster heatmap -------------------------------------------------

lbm_counts <- unique(occ_dt[, .(bldg_id, lbm_cluster)])[, .(lbm_cluster)] |>
  table() |>
  as.data.table()
hourly_occ_lbm <- occ_dt[, .(mean_occ = mean(occ_estimated)),
  by = .(month, day, pseudo_hour, lbm_cluster)
]
hourly_occ_lbm <- merge(hourly_occ_lbm, lbm_counts, by = "lbm_cluster")
hourly_occ_lbm[, lbm_cluster := paste0(lbm_cluster, "\nn = ", N)]

plots[["hourly_occ_lbm"]] <- ggplot() +
  geom_raster(
    data = hourly_occ_lbm,
    aes(x = day, y = pseudo_hour, fill = mean_occ), vjust = 0
  ) +
  scale_x_continuous(breaks = NULL, expand = c(0, 0), name = NULL) +
  scale_y_reverse(name = "hour", expand = c(0, 0), breaks = c(0, 12, 24)) +
  facet_grid(
    cols = vars(month), rows = vars(lbm_cluster),
    scale = "free", space = "free"
  ) +
  scale_fill_gradientn(
    colours  = rev(gradient_occ),
    name     = "Buildings occupied",
    na.value = "transparent",
    limits   = c(0.3, 1),
    oob      = scales::squish,
    breaks   = c(0.3, 0.5, 0.7, 0.9),
    labels   = paste0(c("\U2264 ", rep("", 3)), c(30, 50, 70, 90), "%")
  ) +
  scale_color_manual(values = "black", label = "holidays", name = NULL) +
  guides(
    color = guide_legend(override.aes = list(
      fill = "white",
      key.size = unit(4, "mm")
    )),
    fill = guide_colorbar(title.position = "bottom", title.hjust = 0.5)
  ) +
  theme_time()

plots[["hourly_occ_overview_lbm"]] <- plots[["hourly_occ_lbm"]] +
  plots[["holidays"]] +
  plot_layout(design = shared_layout)

ggsave(
  "plots/03_hourly_occupancy_lbm.pdf",
  plot = plots[["hourly_occ_overview_lbm"]],
  width = 122.5,
  height = 120,
  units = "mm",
  device = pdf
)


## 7.9 Time-cluster calendar plot ---------------------------------------------
#
# Each cell shows the assigned time cluster (T1–T5) for one weekday × week
# combination; holiday periods are overlaid with diagonal hatching.

lbm_date_cls <- occ_dt[bldg_id == bldg_id[1]]
lbm_date_cls[, week := isoweek(time)]
lbm_date_cls[month == "Jan" & week > 10, week := 0] # ISO weeks that cross Jan 1
lbm_date_cls[, holiday := holidays(time)]
lbm_date_cls[, simple_holidays := fifelse(holiday == "No holiday", 0L, 1L)]
lbm_date_cls <- unique(lbm_date_cls[, .(
  week, weekday, lbm_period,
  simple_holidays, month
)])

plots[["lbm_temp"]] <- ggplot(lbm_date_cls, aes(week, weekday)) +
  geom_raster(aes(fill = as.factor(lbm_period))) +
  geom_tile_pattern(
    aes(pattern = as.factor(simple_holidays)),
    fill = NA,
    pattern_density = 0.07,
    pattern_spacing = 0.2,
    pattern_size = 0,
    pattern_fill = "black",
    lineend = "butt"
  ) +
  facet_grid(
    cols = vars(month), rows = vars(weekday),
    scale = "free", space = "free"
  ) +
  scale_fill_manual(values = lbm_temp_colors, name = "Time clusters") +
  scale_pattern_manual(
    name   = NULL,
    values = c("0" = "none", "1" = "stripe"),
    labels = c("1" = "Holidays"),
    breaks = "1"
  ) +
  scale_x_continuous(expand = c(0, 0), name = NULL, breaks = NULL) +
  scale_y_discrete(expand = c(0, 0), name = NULL, breaks = NULL) +
  guides(pattern = guide_legend(
    override.aes = list(pattern_density = 0.2, pattern_spacing = 0.02)
  )) +
  theme_time() +
  theme(panel.spacing.y = unit(0, "mm"), legend.position = "bottom")


## 7.10 Mean occupancy profiles per cluster × time period --------------------

lbm_avg <- occ_dt[, .(mean_occ = mean(occ_estimated)),
  by = .(lbm_cluster, lbm_period, pseudo_hour)
]

plots[["lbm_cluster_mean"]] <- ggplot(
  lbm_avg,
  aes(
    x = pseudo_hour, y = mean_occ,
    group = interaction(lbm_cluster, lbm_period)
  )
) +
  geom_line(linewidth = 0.25) +
  facet_grid(rows = vars(lbm_cluster), cols = vars(lbm_period)) +
  scale_x_continuous(
    breaks = c(0, 6, 12, 18, 24),
    name = "hour of calendar day"
  ) +
  scale_y_continuous(
    name = "avg. occupancy",
    breaks = c(0.4, 0.7, 1.0),
    limits = c(NA, 1.05)
  ) +
  theme_nice() +
  theme(
    text            = element_text(size = 5, colour = "black"),
    axis.text       = element_text(size = 5, colour = "black"),
    panel.spacing.x = unit(1.5, "mm"),
    panel.spacing.y = unit(1, "mm")
  )


## 7.11 Combined FunLBM figure (calendar + profiles) -------------------------

lbm_combined_layout <- c(
  area(t = 1, l = 1, b = 1, r = 6),
  area(t = 2, l = 2, b = 2, r = 5)
)

plots[["lbm_combined"]] <- plots$lbm_temp + plots$lbm_cluster_mean +
  plot_layout(design = lbm_combined_layout) +
  plot_annotation(tag_levels = "a")

ggsave(
  "plots/04_lbm_combined.pdf",
  plot = plots[["lbm_combined"]],
  width = 122.5,
  height = 140,
  units = "mm"
)


# 8. Cross-method cluster comparison ==========================================

## 8.1 Cross-tabulation -------------------------------------------------------
#
# For each pair of CFDA and FunLBM clusters, we compute the proportion of
# buildings shared, normalised once per CFDA cluster and once per FunLBM
# cluster.

cls <- occ_dt[, .SD[1], by = bldg_id]
cls_sum <- cls[, .N, by = .(cfda, lbm_cluster)]
cls_sum[, cfda_norm := N / sum(N), by = cfda]
cls_sum[, lbm_norm := N / sum(N), by = lbm_cluster]
setorder(cls_sum, cfda, lbm_cluster)

bubble_breaks <- c(0.1, seq(0.2, 1, length.out = 5))

# Normalised per CFDA cluster
plots$cls_cfda <- ggplot(
  cls_sum,
  aes(
    x = lbm_cluster, y = cfda,
    size = cfda_norm, color = cfda_norm
  )
) +
  geom_count() +
  scale_color_gradientn(
    labels = percent, colors = rev(gradient_occ),
    na.value = "#FFFFFF", name = "per cfda",
    breaks = bubble_breaks, limit = c(0, 1)
  ) +
  scale_size_continuous(
    labels = percent, name = "per cfda",
    breaks = bubble_breaks, limit = c(0, 1), range = c(2, 9)
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    size = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  coord_fixed() +
  labs(x = "FunLBM", y = "CFDA") +
  theme_nice() +
  theme(legend.title.position = "top", legend.title = element_text(hjust = 0.5))

# Normalised per FunLBM cluster
plots$cls_lbm <- ggplot(
  cls_sum,
  aes(
    x = lbm_cluster, y = cfda,
    size = lbm_norm, color = lbm_norm
  )
) +
  geom_count() +
  scale_color_gradientn(
    labels = percent, colors = rev(gradient_occ),
    na.value = "#FFFFFF", name = "per FunLBM",
    breaks = bubble_breaks, limit = c(0, 1)
  ) +
  scale_size_continuous(
    labels = percent, name = "per FunLBM",
    breaks = bubble_breaks, limit = c(0, 1), range = c(2, 9)
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    size = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  coord_fixed() +
  labs(x = "FunLBM", y = "CFDA") +
  theme_nice() +
  theme(legend.title.position = "top", legend.title = element_text(hjust = 0.5))

ggsave(
  "plots/05_cls_comparison.pdf",
  plot = plots$cls_cfda + plots$cls_lbm + plot_annotation(tag_levels = "a"),
  width = 122.5,
  height = 80,
  units = "mm"
)


## 8.2 Adjusted Rand Index (ARI) ----------------------------------------------
#
# ARI = 1 indicates perfect agreement; ARI ≈ 0 corresponds to chance-level
# agreement.

ari_dt <- occ_dt[, .SD[1], by = bldg_id, .SDcols = c("cfda", "lbm_cluster")]
adjustedRandIndex(ari_dt$cfda, ari_dt$lbm_cluster)
