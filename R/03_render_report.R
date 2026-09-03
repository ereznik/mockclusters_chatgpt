# Render a polished, dependency-light PDF report from analysis_results.rds.
# Run from the project root with:
# Rscript R/03_render_report.R results/analysis_results.rds \
#   output/pdf/mixed_clustering_benchmark.pdf

report_pdf_path <- function(pdf_path) {
  if (length(pdf_path) != 1L || is.na(pdf_path) || !nzchar(pdf_path)) {
    stop("A non-empty PDF filename is required.")
  }
  filename <- basename(pdf_path)
  if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".pdf")
  }
  file.path("output", "pdf", filename)
}

render_benchmark_report <- function(
    results_path = file.path("results", "analysis_results.rds"),
    pdf_path = file.path("output", "pdf", "mixed_clustering_benchmark.pdf")) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The report requires ggplot2. Install it with install.packages('ggplot2').")
  }
  library(grid)
  gg <- asNamespace("ggplot2")
  result <- readRDS(results_path)
  pdf_path <- report_pdf_path(pdf_path)
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  has_truth <- if (!is.null(result$settings$has_truth)) {
    isTRUE(result$settings$has_truth)
  } else {
    !is.null(result$truth) && !all(is.na(result$truth))
  }
  input_label <- if (!is.null(result$settings$input_label)) {
    result$settings$input_label
  } else "input dataset"

  navy <- "#17324D"
  blue <- "#2C7FB8"
  teal <- "#2A9D8F"
  coral <- "#E76F51"
  gold <- "#E9C46A"
  pale <- "#F2F6F8"
  ink <- "#22313F"
  muted <- "#667785"
  method_palette <- c(
    "Gower + PAM" = "#1B9E77",
    "Gower + hierarchical" = "#66A61E",
    "k-prototypes" = "#D95F02",
    "FAMD-style + k-means" = "#7570B3",
    "Latent-embedding mixture" = "#E7298A",
    "Latent class/profile" = "#E6AB02",
    "Gower spectral" = "#A6761D",
    "Grand consensus" = "#17324D"
  )
  base_cluster_colors <- c("#2C7FB8", "#2A9D8F", "#E9C46A", "#E76F51",
                           "#8C6BB1", "#D95F8D", "#66A61E", "#A6761D")
  colors_for_levels <- function(n) {
    if (n <= length(base_cluster_colors)) base_cluster_colors[seq_len(n)]
    else grDevices::hcl.colors(n, "Dark 3")
  }
  cluster_colors <- colors_for_levels(result$settings$k)
  cluster_palette <- stats::setNames(cluster_colors, seq_len(result$settings$k))
  truth_palette <- if (has_truth) {
    truth_levels <- levels(factor(result$truth))
    stats::setNames(colors_for_levels(length(truth_levels)), truth_levels)
  } else NULL
  page_number <- 0L

  theme_report <- function(base_size = 10) {
    gg$theme_minimal(base_size = base_size) +
      gg$theme(
        text = gg$element_text(family = "Helvetica", colour = ink),
        plot.title = gg$element_text(face = "bold", size = base_size + 2,
                                     colour = navy),
        plot.subtitle = gg$element_text(size = base_size - 1, colour = muted),
        panel.grid.minor = gg$element_blank(),
        legend.position = "bottom",
        plot.margin = gg$margin(8, 12, 8, 12)
      )
  }

  footer <- function() {
    grid.lines(x = unit(c(0.055, 0.945), "npc"),
               y = unit(c(0.045, 0.045), "npc"),
               gp = gpar(col = "#D6E0E5", lwd = 0.8))
    grid.text("Mixed-type patient clustering benchmark",
              x = unit(0.055, "npc"), y = unit(0.025, "npc"),
              just = "left", gp = gpar(col = muted, fontsize = 7.5))
    grid.text(as.character(page_number), x = unit(0.945, "npc"),
              y = unit(0.025, "npc"), just = "right",
              gp = gpar(col = muted, fontsize = 7.5))
  }

  begin_page <- function(title, subtitle = NULL) {
    grid.newpage()
    page_number <<- page_number + 1L
    grid.rect(x = 0.5, y = 0.955, width = 1, height = 0.09,
              gp = gpar(fill = navy, col = NA))
    grid.text(title, x = unit(0.055, "npc"), y = unit(0.968, "npc"),
              just = "left", gp = gpar(col = "white", fontsize = 18,
                                        fontface = "bold"))
    if (!is.null(subtitle)) {
      grid.text(subtitle, x = unit(0.055, "npc"), y = unit(0.933, "npc"),
                just = "left", gp = gpar(col = "#DCE8EF", fontsize = 8.5))
    }
    footer()
  }

  wrapped_text <- function(text, x, y, width_chars = 92, fontsize = 10,
                           colour = ink, fontface = "plain", lineheight = 1.25,
                           just = c("left", "top")) {
    lines <- unlist(lapply(text, function(z) strwrap(z, width = width_chars)))
    grid.text(paste(lines, collapse = "\n"), x = unit(x, "npc"),
              y = unit(y, "npc"), just = just,
              gp = gpar(col = colour, fontsize = fontsize,
                        fontface = fontface, lineheight = lineheight))
  }

  bullet_list <- function(items, x = 0.075, y = 0.82, gap = 0.075,
                          width_chars = 88, fontsize = 10) {
    current <- y
    for (item in items) {
      grid.circle(x = unit(x - 0.018, "npc"), y = unit(current - 0.006, "npc"),
                  r = unit(0.004, "npc"), gp = gpar(fill = teal, col = NA))
      lines <- strwrap(item, width = width_chars)
      grid.text(paste(lines, collapse = "\n"), x = unit(x, "npc"),
                y = unit(current, "npc"), just = c("left", "top"),
                gp = gpar(col = ink, fontsize = fontsize, lineheight = 1.2))
      current <- current - gap - max(0, length(lines) - 1L) * 0.025
    }
    invisible(current)
  }

  metric_card <- function(x, y, width, height, label, value, fill = pale,
                          value_colour = navy) {
    grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
                   width = unit(width, "npc"), height = unit(height, "npc"),
                   r = unit(0.012, "npc"), gp = gpar(fill = fill, col = "#D7E2E8"))
    grid.text(value, x = unit(x, "npc"), y = unit(y + 0.012, "npc"),
              gp = gpar(fontsize = 19, fontface = "bold", col = value_colour))
    grid.text(label, x = unit(x, "npc"), y = unit(y - 0.035, "npc"),
              gp = gpar(fontsize = 8, col = muted))
  }

  draw_table <- function(df, top = 0.84, bottom = 0.12,
                         left = 0.055, right = 0.945,
                         font_size = 7.5, header_size = 7.2,
                         widths = NULL) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    nr <- nrow(df) + 1L
    nc <- ncol(df)
    if (is.null(widths)) widths <- rep(1 / nc, nc)
    widths <- widths / sum(widths) * (right - left)
    xs <- left + c(0, cumsum(widths))
    row_h <- (top - bottom) / nr
    for (r in seq_len(nr)) {
      y_top <- top - (r - 1L) * row_h
      fill <- if (r == 1L) navy else if (r %% 2L == 0L) "white" else pale
      grid.rect(x = unit((left + right) / 2, "npc"),
                y = unit(y_top - row_h / 2, "npc"),
                width = unit(right - left, "npc"), height = unit(row_h, "npc"),
                gp = gpar(fill = fill, col = "#D7E2E8", lwd = 0.5))
      for (c in seq_len(nc)) {
        value <- if (r == 1L) names(df)[c] else as.character(df[r - 1L, c])
        grid.text(value,
                  x = unit(xs[c] + 0.007, "npc"),
                  y = unit(y_top - row_h / 2, "npc"),
                  just = "left",
                  gp = gpar(col = if (r == 1L) "white" else ink,
                            fontface = if (r == 1L) "bold" else "plain",
                            fontsize = if (r == 1L) header_size else font_size))
      }
    }
  }

  place_plot <- function(plot, x = 0.06, y = 0.09, width = 0.88, height = 0.80) {
    print(plot, vp = viewport(x = x + width / 2, y = y + height / 2,
                              width = width, height = height))
  }

  draw_raster_matrix <- function(matrix_values, palette, lower, upper,
                                 x, y, width, height) {
    scaled <- (matrix_values - lower) / (upper - lower)
    scaled <- pmax(0, pmin(1, scaled))
    index <- pmax(1L, pmin(length(palette),
                           round(scaled * (length(palette) - 1L)) + 1L))
    colors <- matrix(palette[index], nrow = nrow(matrix_values))
    grid.raster(colors, x = unit(x, "npc"), y = unit(y, "npc"),
                width = unit(width, "npc"), height = unit(height, "npc"),
                interpolate = FALSE)
  }

  draw_cluster_boundaries <- function(labels, left, bottom, width, height,
                                      both_axes = FALSE) {
    sizes <- as.numeric(table(factor(labels, levels = sort(unique(labels)))))
    boundaries <- head(cumsum(sizes) / sum(sizes), -1L)
    for (b in boundaries) {
      y <- bottom + height * (1 - b)
      grid.lines(x = unit(c(left, left + width), "npc"),
                 y = unit(c(y, y), "npc"),
                 gp = gpar(col = coral, lwd = 0.8))
      if (both_axes) {
        x <- left + width * b
        grid.lines(x = unit(c(x, x), "npc"),
                   y = unit(c(bottom, bottom + height), "npc"),
                   gp = gpar(col = coral, lwd = 0.8))
      }
    }
  }

  grDevices::pdf(pdf_path, width = 8.5, height = 11, onefile = TRUE,
                 family = "Helvetica", useDingbats = FALSE, compress = TRUE,
                 title = "Mixed-Type Patient Clustering Benchmark",
                 author = "Reproducible R analysis")
  on.exit(grDevices::dev.off(), add = TRUE)

  # Cover.
  grid.newpage()
  page_number <- page_number + 1L
  grid.rect(gp = gpar(fill = "#F7FAFC", col = NA))
  grid.rect(x = 0.5, y = 0.76, width = 1, height = 0.48,
            gp = gpar(fill = navy, col = NA))
  grid.rect(x = 0.13, y = 0.80, width = 0.012, height = 0.25,
            gp = gpar(fill = teal, col = NA))
  grid.text("Mixed-Type Patient\nClustering Benchmark",
            x = unit(0.17, "npc"), y = unit(0.84, "npc"),
            just = c("left", "top"),
            gp = gpar(col = "white", fontsize = 29, fontface = "bold",
                      lineheight = 1.05))
  grid.text("A CSV-first comparison of eight mixed-data clustering strategies",
            x = unit(0.17, "npc"), y = unit(0.64, "npc"), just = "left",
            gp = gpar(col = "#DCE8EF", fontsize = 11))
  grid.text(sprintf("%d patients  |  %d mixed features  |  %d clusters fitted",
                    result$settings$n, result$settings$p, result$settings$k),
            x = unit(0.17, "npc"), y = unit(0.43, "npc"), just = "left",
            gp = gpar(col = navy, fontsize = 15, fontface = "bold"))
  wrapped_text(
    sprintf("The report analyzes %s using distance-based, prototype-based, dimension-reduced, probabilistic, spectral, and consensus clustering. Dedicated heatmaps show the patient structure produced by every method.", input_label),
    x = 0.17, y = 0.36, width_chars = 70, fontsize = 11, colour = ink
  )
  grid.text(sprintf("Reproducible seed: %s", result$settings$seed),
            x = unit(0.17, "npc"), y = unit(0.17, "npc"), just = "left",
            gp = gpar(col = muted, fontsize = 9))
  footer()

  # Executive summary.
  begin_page("Executive summary", "Headline results from the executed benchmark")
  met <- result$metrics
  cons_idx <- which(met$method == "Grand consensus")
  best_idx <- if (has_truth) which.max(met$ARI) else which.max(met$subsample_stability)
  overall_missing <- mean(result$missingness$missing_rate)
  mean_certainty <- mean(result$membership_certainty)
  metric_card(0.18, 0.82, 0.25, 0.12, "Overall missing cells",
              sprintf("%.1f%%", 100 * overall_missing))
  metric_card(0.50, 0.82, 0.25, 0.12,
              if (has_truth) "Consensus ARI" else "Consensus silhouette",
              sprintf("%.3f", if (has_truth) met$ARI[cons_idx] else met$gower_silhouette[cons_idx]),
              fill = "#E8F5F2", value_colour = teal)
  metric_card(0.82, 0.82, 0.25, 0.12, "Mean membership certainty",
              sprintf("%.1f%%", 100 * mean_certainty), fill = "#FFF6E5", value_colour = "#A86F00")
  first_summary <- if (has_truth) {
    sprintf("The best external recovery was %s (ARI %.3f). The grand consensus achieved ARI %.3f and subsample stability %.3f.",
            met$method[best_idx], met$ARI[best_idx], met$ARI[cons_idx], met$subsample_stability[cons_idx])
  } else {
    sprintf("No truth labels were supplied. %s had the strongest subsample stability (%.3f); the grand consensus stability was %.3f.",
            met$method[best_idx], met$subsample_stability[best_idx], met$subsample_stability[cons_idx])
  }
  last_summary <- if (has_truth) {
    "Known labels are used to score methods only. They are never passed to imputation, dimension reduction, or clustering."
  } else {
    "No external recovery metric is shown because no truth column was supplied. Internal and stability metrics remain available, but external validation is still needed."
  }
  bullet_list(c(
    first_summary,
    sprintf("Across MAR-style imputations, the consensus partitions had mean pairwise ARI %.3f. This isolates uncertainty associated with missing values from ordinary resampling instability.",
            met$imputation_stability[cons_idx]),
    sprintf("Directional MNAR scenarios reassigned between %.1f%% and %.1f%% of patients relative to the baseline scenario.",
            100 * min(result$sensitivity$reassigned_fraction),
            100 * max(result$sensitivity$reassigned_fraction)),
    sprintf("No single internal index proves that %d biological subtypes exist. Retain K=1 as a possibility and require clinical interpretability, resampling stability, and external replication.", result$settings$k),
    last_summary
  ), y = 0.70, gap = 0.085, width_chars = 86, fontsize = 10)

  # Input data overview.
  begin_page("Input data overview", "The feature dictionary controls types, domains, and excluded metadata")
  type_levels <- unique(result$dictionary$type)
  type_summary <- do.call(rbind, lapply(type_levels, function(tp) {
    idx <- result$dictionary$type == tp
    data.frame(
      `Feature type` = tp,
      Features = sum(idx),
      `Mean missing` = sprintf("%.1f%%", 100 * mean(result$missingness$missing_rate[idx])),
      check.names = FALSE
    )
  }))
  draw_table(type_summary, top = 0.84, bottom = 0.55,
             widths = c(0.48, 0.22, 0.30), font_size = 8.5)
  dictionary_note <- if (isTRUE(result$settings$dictionary_inferred)) {
    "Feature types were inferred from the CSV and saved to feature_dictionary_used.csv. Review ordinal variables and clinical domains before treating this as a final analysis."
  } else {
    "A supplied feature dictionary determined variable types and domain weights. Columns absent from the dictionary were not used for clustering."
  }
  bullet_list(c(
    sprintf("Input: %s. The analysis used %d rows and %d clustering features.", input_label, result$settings$n, result$settings$p),
    dictionary_note,
    sprintf("The dictionary contains %d weighting domains. Each domain receives equal total Gower weight, so large feature blocks do not dominate solely by column count.", length(unique(result$dictionary$domain))),
    sprintf("Observed missingness is %.1f%% overall. Stochastic mixed-type hot-deck imputations preserve observed categories and propagate assignment uncertainty.", 100 * overall_missing),
    if (has_truth) "A truth column was supplied for external evaluation and excluded from all fitting steps." else "No truth column was supplied; ARI, NMI, and label-matched accuracy are therefore intentionally unavailable."
  ), y = 0.47, gap = 0.060, width_chars = 94, fontsize = 9.1)

  # Missingness.
  begin_page("Observed missingness", "Feature-level rates; colors indicate declared measurement type")
  miss <- result$missingness
  miss <- miss[order(miss$missing_rate, decreasing = TRUE), ]
  miss_top <- head(miss, 30L)
  miss_top$feature <- factor(miss_top$feature, levels = rev(miss_top$feature))
  p_miss <- gg$ggplot(miss_top, gg$aes(x = missing_rate, y = feature, fill = type)) +
    gg$geom_col(width = 0.72) +
    gg$scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                          expand = gg$expansion(mult = c(0, 0.06))) +
    gg$labs(x = "Missing rate", y = NULL,
            title = "Thirty features with the most missing observations",
            subtitle = sprintf("Overall rate: %.1f%%; maximum feature rate: %.1f%%",
                               100 * mean(miss$missing_rate),
                               100 * max(miss$missing_rate)), fill = "Type") +
    theme_report(9) +
    gg$theme(axis.text.y = gg$element_text(size = 7.2))
  place_plot(p_miss, y = 0.08, height = 0.82)

  # Methods table.
  begin_page("Approaches implemented", "All methods use the same type dictionary and imputed datasets")
  approach_table <- data.frame(
    Approach = c("Gower + PAM", "Gower + hierarchical", "k-prototypes",
                 "FAMD-style + k-means", "Latent-embedding mixture",
                 "Latent class/profile", "Gower spectral", "Grand consensus"),
    Representation = c("Weighted mixed dissimilarity", "Weighted mixed dissimilarity",
                       "Standardized numeric + category modes", "Balanced mixed latent axes",
                       "Balanced mixed latent axes", "Typed conditional distributions",
                       "Gower affinity eigenvectors", "Co-assignment matrix"),
    Assignment = c("Medoid", "Average-linkage cut", "Prototype", "Centroid",
                   "Posterior mode", "Posterior mode", "Spectral centroid", "Consensus medoid"),
    Key_caveat = c("Feature weighting", "Linkage sensitivity", "Lambda and local minima",
                   "Two-stage approximation", "Latent Gaussian assumption",
                   "Conditional independence", "Kernel bandwidth", "Can average weak answers"),
    check.names = FALSE
  )
  names(approach_table)[4] <- "Key caveat"
  draw_table(approach_table, top = 0.86, bottom = 0.19,
             widths = c(0.24, 0.30, 0.19, 0.27), font_size = 7.7)
  wrapped_text(
    "Implementation note: the FAMD-style embedding is a dependency-light PCAmix approximation. It standardizes numeric/ordinal variables and balances each nominal indicator block before PCA. Production work can replace it with FactoMineR::FAMD or PCAmixdata without changing the evaluation harness.",
    x = 0.06, y = 0.14, width_chars = 105, fontsize = 8.4, colour = muted
  )

  # Performance table.
  begin_page("Method performance", "Consensus per method across eight imputations; higher is better except runtime")
  perf <- result$metrics
  perf_display <- data.frame(
    Method = perf$method,
    ARI = sprintf("%.3f", perf$ARI),
    NMI = sprintf("%.3f", perf$NMI),
    Accuracy = sprintf("%.1f%%", 100 * perf$best_match_accuracy),
    Silhouette = sprintf("%.3f", perf$gower_silhouette),
    `MI stability` = sprintf("%.3f", perf$imputation_stability),
    `Subsample stability` = sprintf("%.3f", perf$subsample_stability),
    `Runtime (s)` = sprintf("%.2f", perf$mean_runtime_seconds),
    check.names = FALSE
  )
  draw_table(perf_display, top = 0.85, bottom = 0.28,
             widths = c(0.27, 0.075, 0.075, 0.10, 0.11, 0.12, 0.15, 0.10),
             font_size = 6.9, header_size = 6.5)
  metric_explanation <- if (has_truth) {
    "ARI and NMI compare clusters to the supplied truth labels. Accuracy is calculated after the best one-to-one label permutation."
  } else {
    "ARI, NMI, and accuracy are NA because no truth labels were supplied. This is expected for ordinary unsupervised analysis."
  }
  wrapped_text(
    paste(metric_explanation, "Silhouette uses the first completed-data Gower distance. MI stability is mean pairwise ARI across imputation-specific partitions. Subsample stability compares refitted 80% subsets with the full-data solution on the same patients."),
    x = 0.06, y = 0.22, width_chars = 108, fontsize = 8.7, colour = muted
  )
  wrapped_text(
    "Runtime values are indicative on this machine and omit report rendering. Shared preprocessing makes exact cross-method timing comparisons approximate.",
    x = 0.06, y = 0.12, width_chars = 108, fontsize = 8.7, colour = muted
  )

  # Performance plot.
  begin_page(if (has_truth) "Recovery of supplied labels" else "Internal performance and stability",
             if (has_truth) "External metrics are used only for evaluation" else "No external truth labels were supplied")
  metric_long <- if (has_truth) {
    rbind(
      data.frame(method = perf$method, metric = "ARI", value = perf$ARI),
      data.frame(method = perf$method, metric = "NMI", value = perf$NMI),
      data.frame(method = perf$method, metric = "Best-match accuracy", value = perf$best_match_accuracy)
    )
  } else {
    rbind(
      data.frame(method = perf$method, metric = "Gower silhouette", value = perf$gower_silhouette),
      data.frame(method = perf$method, metric = "MI stability", value = perf$imputation_stability),
      data.frame(method = perf$method, metric = "Subsample stability", value = perf$subsample_stability)
    )
  }
  metric_long$method <- factor(metric_long$method, levels = rev(perf$method))
  plot_colors <- if (has_truth) {
    c("ARI" = blue, "NMI" = teal, "Best-match accuracy" = coral)
  } else {
    c("Gower silhouette" = blue, "MI stability" = teal, "Subsample stability" = coral)
  }
  p_perf <- gg$ggplot(metric_long,
                      gg$aes(x = value, y = method, fill = metric)) +
    gg$geom_col(position = gg$position_dodge(width = 0.72), width = 0.66) +
    gg$geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed",
                  colour = "#B8C4CB", linewidth = 0.4) +
    gg$scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    gg$scale_fill_manual(values = plot_colors) +
    gg$labs(x = "Score", y = NULL, fill = NULL,
            title = if (has_truth) "Agreement with supplied labels" else "Internal separation and reproducibility") +
    theme_report(10)
  place_plot(p_perf, y = 0.11, height = 0.76)

  # Stability plot.
  begin_page("Stability and missing-data uncertainty", "Strong recovery is useful only when assignments survive perturbation")
  perf$point_size <- if (has_truth) perf$ARI else perf$gower_silhouette
  size_title <- if (has_truth) "ARI vs truth" else "Gower silhouette"
  size_limits <- if (has_truth) c(0, 1) else range(c(0, perf$point_size), na.rm = TRUE)
  if (diff(size_limits) < 1e-8) size_limits <- c(0, max(0.1, size_limits[2L]))
  p_stab <- gg$ggplot(perf,
                      gg$aes(x = imputation_stability, y = subsample_stability,
                             size = point_size, colour = method)) +
    gg$geom_hline(yintercept = 0.80, linetype = "dashed", colour = "#C4CDD2") +
    gg$geom_vline(xintercept = 0.80, linetype = "dashed", colour = "#C4CDD2") +
    gg$geom_point(alpha = 0.90) +
    gg$scale_colour_manual(values = method_palette) +
    gg$scale_size_continuous(range = c(3.5, 8), limits = size_limits) +
    gg$coord_cartesian(xlim = c(min(0.35, min(perf$imputation_stability) - 0.03), 1.01),
                       ylim = c(min(0.35, min(perf$subsample_stability) - 0.03), 1.01)) +
    gg$labs(x = "Stability across imputations (pairwise ARI)",
            y = "Repeated 80% subsample stability (ARI)",
            colour = NULL, size = size_title,
            title = "Two distinct stability questions") +
    theme_report(10) +
    gg$theme(legend.position = "right",
             legend.text = gg$element_text(size = 6.3),
             legend.title = gg$element_text(size = 7),
             legend.key.height = unit(0.030, "npc")) +
    gg$guides(colour = gg$guide_legend(order = 1, override.aes = list(size = 4)),
              size = gg$guide_legend(order = 2))
  place_plot(p_stab, y = 0.11, height = 0.77)

  # K selection.
  k_selection_subtitle <- if (identical(result$settings$k_selection_mode, "auto")) {
    sprintf("Automatic rank aggregation selected K=%d; no single criterion should decide K",
            result$settings$k)
  } else {
    sprintf("K=%d was supplied; diagnostics across the candidate range remain essential",
            result$settings$k)
  }
  begin_page("Choosing the number of clusters", k_selection_subtitle)
  ks <- result$k_selection
  p_k1 <- gg$ggplot(ks, gg$aes(k, gower_pam_silhouette)) +
    gg$geom_line(colour = blue, linewidth = 0.8) + gg$geom_point(colour = blue, size = 2) +
    gg$geom_vline(xintercept = result$settings$k, linetype = "dashed", colour = coral) +
    gg$scale_x_continuous(breaks = ks$k) +
    gg$labs(x = "K", y = "Silhouette", title = "Gower + PAM") + theme_report(8)
  p_k2 <- gg$ggplot(ks, gg$aes(k, famd_silhouette)) +
    gg$geom_line(colour = teal, linewidth = 0.8) + gg$geom_point(colour = teal, size = 2) +
    gg$geom_vline(xintercept = result$settings$k, linetype = "dashed", colour = coral) +
    gg$scale_x_continuous(breaks = ks$k) +
    gg$labs(x = "K", y = "Silhouette", title = "FAMD-style + k-means") + theme_report(8)
  p_k3 <- gg$ggplot(ks, gg$aes(k, embedding_gmm_bic)) +
    gg$geom_line(colour = "#8C6BB1", linewidth = 0.8) + gg$geom_point(colour = "#8C6BB1", size = 2) +
    gg$geom_vline(xintercept = result$settings$k, linetype = "dashed", colour = coral) +
    gg$scale_x_continuous(breaks = ks$k) +
    gg$labs(x = "K", y = "BIC (lower is better)", title = "Latent mixture") + theme_report(8)
  print(p_k1, vp = viewport(x = 0.5, y = 0.74, width = 0.88, height = 0.25))
  print(p_k2, vp = viewport(x = 0.5, y = 0.47, width = 0.88, height = 0.25))
  print(p_k3, vp = viewport(x = 0.5, y = 0.19, width = 0.88, height = 0.25))

  # Embedding.
  begin_page("Mixed-data latent view", "First two FAMD-style axes provide a common visual reference")
  left_labels <- if (has_truth) factor(result$truth) else factor(result$method_consensus[, "Gower + PAM"])
  emb_df <- data.frame(
    Axis1 = result$embedding$scores[, 1L],
    Axis2 = result$embedding$scores[, 2L],
    Left = left_labels,
    Consensus = factor(result$grand_consensus)
  )
  p_e1 <- gg$ggplot(emb_df, gg$aes(Axis1, Axis2, colour = Left)) +
    gg$geom_point(alpha = 0.70, size = 1.5) +
    gg$scale_colour_manual(values = if (has_truth) truth_palette else cluster_palette) +
    gg$labs(title = if (has_truth) "Colored by supplied truth" else "Colored by Gower + PAM",
            x = "Axis 1", y = "Axis 2", colour = if (has_truth) "Truth" else "Cluster") +
    theme_report(8)
  p_e2 <- gg$ggplot(emb_df, gg$aes(Axis1, Axis2, colour = Consensus)) +
    gg$geom_point(alpha = 0.70, size = 1.5) +
    gg$scale_colour_manual(values = unname(cluster_palette)) +
    gg$labs(title = "Colored by grand consensus", x = "Axis 1", y = "Axis 2", colour = "Cluster") +
    theme_report(8)
  print(p_e1, vp = viewport(x = 0.27, y = 0.49, width = 0.44, height = 0.76))
  print(p_e2, vp = viewport(x = 0.73, y = 0.49, width = 0.44, height = 0.76))

  # Co-assignment heatmap.
  begin_page("Consensus co-assignment", "Cell color is the fraction of method-imputation partitions placing two patients together")
  ord <- if (has_truth) {
    order(result$grand_consensus, result$truth)
  } else {
    order(result$grand_consensus, result$embedding$scores[, 1L])
  }
  co <- result$coassignment[ord, ord]
  pal <- grDevices::colorRampPalette(c("#F7FAFC", "#BFDCEB", "#2C7FB8", "#17324D"))(256)
  image_index <- pmax(1L, pmin(256L, round(co * 255L) + 1L))
  grid.raster(matrix(pal[image_index], nrow = nrow(co)),
              x = unit(0.50, "npc"), y = unit(0.50, "npc"),
              width = unit(0.82, "npc"), height = unit(0.72, "npc"),
              interpolate = FALSE)
  sizes <- as.numeric(table(result$grand_consensus))
  bounds <- cumsum(sizes) / sum(sizes)
  for (b in head(bounds, -1L)) {
    x <- 0.09 + 0.82 * b
    y <- 0.14 + 0.72 * (1 - b)
    grid.lines(x = unit(c(x, x), "npc"), y = unit(c(0.14, 0.86), "npc"),
               gp = gpar(col = coral, lwd = 1.0))
    grid.lines(x = unit(c(0.09, 0.91), "npc"), y = unit(c(y, y), "npc"),
               gp = gpar(col = coral, lwd = 1.0))
  }
  grid.text(if (has_truth) "Patients sorted by consensus cluster, then supplied truth" else "Patients sorted by consensus cluster, then latent axis 1",
            x = unit(0.50, "npc"), y = unit(0.10, "npc"),
            gp = gpar(col = muted, fontsize = 8.5))

  # One dedicated dual-heatmap page for every method and the grand consensus.
  heat_methods <- c(result$method_names, "Grand consensus")
  design <- result$embedding$design
  if (is.null(design)) stop("Heatmap-ready design matrix is absent from the analysis result.")
  co_palette <- grDevices::colorRampPalette(c("#F7FAFC", "#BFDCEB", "#2C7FB8", "#17324D"))(256)
  feature_palette <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(256)
  for (method in heat_methods) {
    labels <- if (method == "Grand consensus") {
      result$grand_consensus
    } else result$method_consensus[, method]
    co_method <- if (method == "Grand consensus") {
      result$coassignment
    } else result$method_coassignment[[method]]
    row_order <- order(labels, result$embedding$scores[, 1L])
    ordered_labels <- labels[row_order]

    effect <- vapply(seq_len(ncol(design)), function(j) {
      x <- design[, j]
      total <- sum((x - mean(x))^2)
      if (!is.finite(total) || total < 1e-10) return(0)
      group_means <- tapply(x, labels, mean)
      group_sizes <- table(labels)
      sum(as.numeric(group_sizes[names(group_means)]) *
            (group_means - mean(x))^2) / total
    }, numeric(1L))
    top_columns <- head(order(effect, decreasing = TRUE), min(32L, ncol(design)))
    feature_heat <- design[row_order, top_columns, drop = FALSE]
    feature_heat[!is.finite(feature_heat)] <- 0
    co_heat <- co_method[row_order, row_order]
    metric_row <- result$metrics[result$metrics$method == method, , drop = FALSE]

    begin_page(paste("Method heatmaps:", method),
               "Left: stability across imputations. Right: patients by the most cluster-associated encoded features")
    metric_card(0.16, 0.83, 0.20, 0.10,
                if (has_truth) "ARI" else "Silhouette",
                sprintf("%.3f", if (has_truth) metric_row$ARI else metric_row$gower_silhouette),
                fill = "#E8F5F2", value_colour = teal)
    metric_card(0.39, 0.83, 0.20, 0.10, "MI stability",
                sprintf("%.3f", metric_row$imputation_stability))
    metric_card(0.62, 0.83, 0.20, 0.10, "Subsample stability",
                sprintf("%.3f", metric_row$subsample_stability),
                fill = "#FFF6E5", value_colour = "#A86F00")
    metric_card(0.85, 0.83, 0.20, 0.10, "Gower silhouette",
                sprintf("%.3f", metric_row$gower_silhouette))

    grid.text("Patient-patient co-assignment", x = unit(0.27, "npc"),
              y = unit(0.735, "npc"), gp = gpar(fontsize = 10, fontface = "bold", col = navy))
    grid.text(sprintf("Patient-feature heatmap (%d encoded columns)", length(top_columns)),
              x = unit(0.73, "npc"), y = unit(0.735, "npc"),
              gp = gpar(fontsize = 10, fontface = "bold", col = navy))

    draw_raster_matrix(co_heat, co_palette, 0, 1,
                       x = 0.27, y = 0.43, width = 0.40, height = 0.54)
    draw_raster_matrix(feature_heat, feature_palette, -2.5, 2.5,
                       x = 0.73, y = 0.43, width = 0.40, height = 0.54)
    draw_cluster_boundaries(ordered_labels, left = 0.07, bottom = 0.16,
                            width = 0.40, height = 0.54, both_axes = TRUE)
    draw_cluster_boundaries(ordered_labels, left = 0.53, bottom = 0.16,
                            width = 0.40, height = 0.54, both_axes = FALSE)
    grid.rect(x = unit(0.27, "npc"), y = unit(0.43, "npc"),
              width = unit(0.40, "npc"), height = unit(0.54, "npc"),
              gp = gpar(fill = NA, col = "#AFC0C9", lwd = 0.7))
    grid.rect(x = unit(0.73, "npc"), y = unit(0.43, "npc"),
              width = unit(0.40, "npc"), height = unit(0.54, "npc"),
              gp = gpar(fill = NA, col = "#AFC0C9", lwd = 0.7))

    sizes_text <- paste(sprintf("C%s: %d", names(table(labels)), as.numeric(table(labels))),
                        collapse = "   |   ")
    grid.text(paste("Cluster sizes:", sizes_text), x = unit(0.50, "npc"),
              y = unit(0.115, "npc"), gp = gpar(fontsize = 8.2, col = ink))
    top_names <- colnames(design)[top_columns]
    top_names <- gsub("::", "=", top_names, fixed = TRUE)
    wrapped_text(paste("Most associated encoded columns:", paste(head(top_names, 8L), collapse = ", ")),
                 x = 0.07, y = 0.085, width_chars = 115, fontsize = 7.2, colour = muted)
  }

  # Profiles.
  begin_page("Consensus cluster profiles", "Mean standardized values for selected numeric or ordinal features")
  prof_long <- as.data.frame(as.table(result$profile_means), stringsAsFactors = FALSE)
  names(prof_long) <- c("Cluster", "Feature", "Mean")
  prof_long$Feature <- factor(prof_long$Feature,
                              levels = rev(colnames(result$profile_means)))
  p_prof <- gg$ggplot(prof_long, gg$aes(x = Cluster, y = Feature, fill = Mean)) +
    gg$geom_tile(colour = "white", linewidth = 0.3) +
    gg$scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D95F02",
                            midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish) +
    gg$labs(x = NULL, y = NULL, fill = "Z mean",
            title = "Profiles describe the fitted groups; they are not independent validation") +
    theme_report(9) +
    gg$theme(axis.text.y = gg$element_text(size = 7.3))
  place_plot(p_prof, y = 0.10, height = 0.78)

  # Certainty and confusion.
  begin_page("Membership uncertainty", "Consensus support exposes patients whose assignments are not robust")
  certainty_df <- data.frame(certainty = result$membership_certainty,
                             cluster = factor(result$grand_consensus))
  p_cert <- gg$ggplot(certainty_df, gg$aes(certainty, fill = cluster)) +
    gg$geom_histogram(binwidth = 0.04, boundary = 0, colour = "white") +
    gg$scale_fill_manual(values = cluster_palette) +
    gg$scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")) +
    gg$labs(x = "Fraction of aligned partitions supporting consensus label",
            y = "Patients", fill = "Cluster", title = "Distribution of membership certainty") +
    theme_report(8)
  print(p_cert, vp = viewport(x = 0.5, y = 0.69, width = 0.88, height = 0.37))
  low <- order(result$membership_certainty)[seq_len(10L)]
  uncertain <- if (has_truth) {
    data.frame(
      Patient = result$patient_id[low],
      Truth = as.character(result$truth[low]),
      Consensus = result$grand_consensus[low],
      Certainty = sprintf("%.1f%%", 100 * result$membership_certainty[low]),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Patient = result$patient_id[low],
      Consensus = result$grand_consensus[low],
      Certainty = sprintf("%.1f%%", 100 * result$membership_certainty[low]),
      stringsAsFactors = FALSE
    )
  }
  grid.text("Ten least-certain patients", x = unit(0.06, "npc"), y = unit(0.46, "npc"),
            just = "left", gp = gpar(fontsize = 10, fontface = "bold", col = navy))
  draw_table(uncertain, top = 0.43, bottom = 0.10,
             widths = if (has_truth) c(0.30, 0.22, 0.25, 0.23) else c(0.45, 0.28, 0.27),
             font_size = 7.6)

  # MNAR sensitivity.
  mnar_count <- if (!is.null(result$settings$mnar_feature_count)) result$settings$mnar_feature_count else 0L
  begin_page("MNAR sensitivity analysis",
             sprintf("Only originally missing values in %d dictionary-designated features are shifted", mnar_count))
  sens <- result$sensitivity
  sens$Scenario <- factor(sens$interpretation, levels = sens$interpretation)
  sens_long <- if (has_truth) {
    rbind(
      data.frame(Scenario = sens$Scenario, Metric = "ARI vs truth", Value = sens$ARI_vs_truth),
      data.frame(Scenario = sens$Scenario, Metric = "ARI vs baseline", Value = sens$ARI_vs_baseline)
    )
  } else {
    data.frame(Scenario = sens$Scenario, Metric = "ARI vs baseline", Value = sens$ARI_vs_baseline)
  }
  p_sens <- gg$ggplot(sens_long, gg$aes(Scenario, Value, fill = Metric)) +
    gg$geom_col(position = gg$position_dodge(width = 0.72), width = 0.65) +
    gg$scale_fill_manual(values = c("ARI vs truth" = teal, "ARI vs baseline" = navy)) +
    gg$scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    gg$labs(x = NULL, y = "Adjusted Rand index", fill = NULL,
            title = "Sensitivity of pooled representative methods") +
    theme_report(8) +
    gg$theme(axis.text.x = gg$element_text(angle = 15, hjust = 1))
  print(p_sens, vp = viewport(x = 0.5, y = 0.67, width = 0.88, height = 0.42))
  sens_display <- if (has_truth) {
    data.frame(
      Delta = sprintf("%+.0f", sens$delta),
      Scenario = sens$interpretation,
      `ARI vs truth` = sprintf("%.3f", sens$ARI_vs_truth),
      `ARI vs baseline` = sprintf("%.3f", sens$ARI_vs_baseline),
      Reassigned = sprintf("%.1f%%", 100 * sens$reassigned_fraction),
      check.names = FALSE
    )
  } else {
    data.frame(
      Delta = sprintf("%+.0f", sens$delta),
      Scenario = sens$interpretation,
      `ARI vs baseline` = sprintf("%.3f", sens$ARI_vs_baseline),
      Reassigned = sprintf("%.1f%%", 100 * sens$reassigned_fraction),
      check.names = FALSE
    )
  }
  draw_table(sens_display, top = 0.43, bottom = 0.24,
             widths = if (has_truth) c(0.10, 0.35, 0.18, 0.20, 0.17) else c(0.12, 0.43, 0.25, 0.20),
             font_size = 7.4)
  wrapped_text(
    if (mnar_count > 0L) {
      "Delta is measured in observed standard deviations for numeric features and approximately one category step for ordinal features. These scenarios illustrate sensitivity, not identification: the true MNAR mechanism cannot be recovered from observed data alone."
    } else {
      "No mnar_direction values were specified, so the three scenarios are identical. Edit feature_dictionary_used.csv to mark plausible directions as -1 or +1, then rerun the CSV command to obtain an informative sensitivity analysis."
    },
    x = 0.06, y = 0.18, width_chars = 108, fontsize = 8.6, colour = muted
  )

  # Interpretation and limitations.
  begin_page("How to interpret this benchmark", "What transfers to a real cohort and what does not")
  bullet_list(c(
    "External recovery scores transfer only as evidence that the implementations work on this simulation. A real cohort has no hidden truth labels, and good silhouette values do not establish biological subtypes.",
    "The stochastic nearest-neighbor hot-deck procedure preserves mixed-data values and propagates some imputation uncertainty. For publication-grade analysis, consider a carefully specified MICE or joint-model imputation and repeat the same consensus workflow.",
    "The native latent class/profile model uses diagonal Gaussian distributions for continuous/log-count variables and multinomials for ordinal and nominal variables. It is intentionally parsimonious; specialized count and ordinal likelihoods may be preferable.",
    "The FAMD-style implementation is transparent and dependency-light but is not a drop-in reproduction of every weighting convention used by FactoMineR or PCAmixdata.",
    "Cluster-feature differences are descriptive because those features constructed the clusters. Validate with prespecified variables not used for clustering, then replicate in an external cohort.",
    "If the scientific goal is prediction, treatment selection, or outcome risk, a supervised model may be more appropriate than clustering. Clusters should not be interpreted causally."
  ), y = 0.84, gap = 0.076, width_chars = 94, fontsize = 9.1)

  # Rerun page.
  begin_page("Rerunning with your own data", "The project is designed to be edited rather than treated as a black box")
  wrapped_text("Fastest complete rerun", x = 0.06, y = 0.86,
               fontsize = 12, fontface = "bold", colour = navy)
  grid.roundrect(x = 0.5, y = 0.80, width = 0.88, height = 0.085,
                 gp = gpar(fill = "#EEF3F6", col = "#D7E2E8"))
  grid.text("Rscript analyze_csv.R --input path/to/patients.csv --k auto",
            x = unit(0.08, "npc"), y = unit(0.80, "npc"),
            just = "left", gp = gpar(fontfamily = "Courier", fontsize = 11, col = ink))
  bullet_list(c(
    "Put one patient per row. Use --id-column for an identifier, --truth-column only for optional external labels, and --exclude for outcomes, site, batch, or other columns that should not define clusters.",
    "A dictionary is inferred when --dictionary is omitted and saved as feature_dictionary_used.csv. Review it, especially numeric ordinal variables, then rerun with --dictionary pointing to the edited file.",
    "Allowed types are continuous, discrete, ordinal, categorical, and binary. The domain column controls block weighting. Optional role and mnar_direction columns support reporting and sensitivity analysis.",
    "Use --k auto for a rank aggregation of two silhouette diagnostics and latent-mixture BIC, or provide an integer. Always inspect K-selection, stability, cluster sizes, uncertainty, and clinical coherence together.",
    "Increase m_imputations and n_subsamples for a final analysis. Eight and ten are intentionally modest defaults for an executable demonstration.",
    "Machine-readable outputs include method metrics, cluster assignments, membership certainty, missingness, K selection, MNAR sensitivity, resampling stability, the dictionary used, and a run manifest."
  ), y = 0.71, gap = 0.075, width_chars = 92, fontsize = 9.1)

  # Output inventory.
  begin_page("Project output inventory", "Files produced by the reproducible pipeline")
  inventory <- data.frame(
    File = c("feature_dictionary_used.csv", "auto_k_diagnostics.csv (if --k auto)",
             "run_manifest.csv", "results/method_metrics.csv",
             "results/cluster_assignments.csv", "results/membership_certainty.csv",
             "results/k_selection.csv", "results/mnar_sensitivity.csv",
             "results/analysis_results.rds", "output/pdf/*_clustering_report.pdf"),
    Purpose = c("Feature types and domains actually used", "Diagnostics used when --k auto is selected",
                "Input, output, K, ID, and truth settings", "External, internal, stability, and runtime metrics",
                "Consensus and method-specific patient labels", "Patient-level assignment support",
                "Candidate-K diagnostics", "Directional departure-from-MAR scenarios",
                "Complete report-ready analysis object", "This rendered benchmark report"),
    stringsAsFactors = FALSE
  )
  draw_table(inventory, top = 0.86, bottom = 0.18,
             widths = c(0.43, 0.57), font_size = 8.1)
  wrapped_text(
    "Machine-readable files are written below the chosen --output-dir. PDF reports are always written to output/pdf; --pdf changes only the filename.",
    x = 0.06, y = 0.13, width_chars = 105, fontsize = 8.8, colour = muted
  )

  # Reproducibility.
  begin_page("Reproducibility record", "Software information captured during execution")
  session_lines <- result$session_info
  if (length(session_lines) > 46L) session_lines <- session_lines[seq_len(46L)]
  grid.text(paste(session_lines, collapse = "\n"),
            x = unit(0.06, "npc"), y = unit(0.87, "npc"),
            just = c("left", "top"),
            gp = gpar(fontfamily = "Courier", fontsize = 7.1,
                      col = ink, lineheight = 1.12))
  grid.text(sprintf("Analysis timestamp (UTC): %s", result$settings$generated_at),
            x = unit(0.06, "npc"), y = unit(0.10, "npc"), just = "left",
            gp = gpar(col = muted, fontsize = 8.5))

  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  message("PDF report written to ", pdf_path)
  invisible(pdf_path)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  results_path <- if (length(args) >= 1L) args[1L] else file.path("results", "analysis_results.rds")
  pdf_path <- if (length(args) >= 2L) args[2L] else file.path("output", "pdf", "mixed_clustering_benchmark.pdf")
  render_benchmark_report(results_path, pdf_path)
}
