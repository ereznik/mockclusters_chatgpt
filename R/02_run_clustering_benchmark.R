# Run all mixed-data clustering approaches and save machine-readable results.
# Standalone use:
# Rscript R/02_run_clustering_benchmark.R data/mock_patients_typed.rds \
#   data/feature_dictionary.csv

source(file.path("R", "clustering_functions.R"))

run_clustering_benchmark <- function(
    data_path = file.path("data", "mock_patients_typed.rds"),
    dictionary_path = file.path("data", "feature_dictionary.csv"),
    results_dir = "results",
    seed = 20260902L,
    k = 4L,
    m_imputations = 8L,
    n_subsamples = 10L,
    id_column = NULL,
    truth_column = NULL) {

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dictionary <- utils::read.csv(dictionary_path, stringsAsFactors = FALSE)
  required_dict <- c("feature", "type", "domain")
  if (!all(required_dict %in% names(dictionary))) {
    stop("Dictionary must contain: ", paste(required_dict, collapse = ", "))
  }
  if (!"mnar_direction" %in% names(dictionary)) dictionary$mnar_direction <- 0L
  if (!"role" %in% names(dictionary)) dictionary$role <- "unspecified"

  if (grepl("\\.rds$", data_path, ignore.case = TRUE)) {
    data <- readRDS(data_path)
  } else {
    data <- utils::read.csv(data_path, stringsAsFactors = FALSE,
                            na.strings = c("", "NA"), check.names = FALSE)
  }
  data <- restore_types(data, dictionary)
  n <- nrow(data)
  if (is.null(id_column) && "patient_id" %in% names(data)) id_column <- "patient_id"
  truth <- if (!is.null(truth_column) && truth_column %in% names(data)) {
    data[[truth_column]]
  } else rep(NA, n)
  patient_id <- if (!is.null(id_column) && id_column %in% names(data)) {
    as.character(data[[id_column]])
  } else sprintf("row_%04d", seq_len(n))
  if (anyNA(patient_id) || anyDuplicated(patient_id)) {
    stop("The ID column must be complete and unique.")
  }

  message("Creating stochastic mixed-type imputations...")
  imputations <- stochastic_hotdeck_impute(
    data, dictionary, m = m_imputations, seed = seed + 100L
  )

  method_names <- c("Gower + PAM", "Gower + hierarchical", "k-prototypes",
                    "FAMD-style + k-means", "Latent-embedding mixture",
                    "Latent class/profile", "Gower spectral")
  partition_array <- array(
    NA_integer_,
    dim = c(n, length(method_names), m_imputations),
    dimnames = list(patient_id, method_names, paste0("imp_", seq_len(m_imputations)))
  )
  runtime_array <- matrix(
    NA_real_, nrow = length(method_names), ncol = m_imputations,
    dimnames = list(method_names, paste0("imp_", seq_len(m_imputations)))
  )
  first_details <- NULL

  for (imp in seq_len(m_imputations)) {
    message(sprintf("  Imputation %d/%d: running seven base approaches...",
                    imp, m_imputations))
    fit <- run_all_methods(
      imputations[[imp]], dictionary, k = k,
      seed = seed + 1000L * imp, keep_details = imp == 1L
    )
    for (method in method_names) {
      partition_array[, method, imp] <- fit$labels[[method]]
      runtime_array[method, imp] <- fit$runtimes[[method]]
    }
    if (imp == 1L) first_details <- fit$details
  }

  # Pool each method across imputations, then pool every method/imputation together.
  method_consensus <- matrix(
    NA_integer_, n, length(method_names),
    dimnames = list(patient_id, method_names)
  )
  imputation_stability <- setNames(numeric(length(method_names)), method_names)
  method_coassign <- vector("list", length(method_names))
  names(method_coassign) <- method_names
  for (method in method_names) {
    parts <- partition_array[, method, , drop = FALSE]
    dim(parts) <- c(n, m_imputations)
    cons <- consensus_from_partitions(parts, k = k)
    method_consensus[, method] <- cons$labels
    method_coassign[[method]] <- cons$coassignment
    imputation_stability[method] <- pairwise_partition_ari(parts)
  }
  all_partitions <- matrix(partition_array, nrow = n)
  grand <- consensus_from_partitions(all_partitions, k = k)

  per_imputation_consensus <- matrix(NA_integer_, n, m_imputations)
  for (imp in seq_len(m_imputations)) {
    per_imputation_consensus[, imp] <- consensus_from_partitions(
      partition_array[, , imp], k = k
    )$labels
  }
  grand_imputation_stability <- pairwise_partition_ari(per_imputation_consensus)

  # Patient-level certainty: fraction of aligned base partitions supporting the
  # grand-consensus label.
  aligned_partitions <- apply(all_partitions, 2L, align_labels,
                              reference = grand$labels)
  grand_code <- as.integer(factor(grand$labels))
  membership_certainty <- rowMeans(aligned_partitions == grand_code)

  # Subsampling stability repeats every algorithm on 80% of patients.
  message(sprintf("Estimating stability with %d repeated 80%% subsamples...",
                  n_subsamples))
  stability_values <- matrix(
    NA_real_, nrow = length(method_names) + 1L, ncol = n_subsamples,
    dimnames = list(c(method_names, "Grand consensus"), paste0("sub_", seq_len(n_subsamples)))
  )
  set.seed(seed + 200L)
  for (b in seq_len(n_subsamples)) {
    idx <- sort(sample.int(n, floor(0.80 * n), replace = FALSE))
    fit_b <- run_all_methods(
      imputations[[1L]][idx, , drop = FALSE], dictionary, k = k,
      seed = seed + 5000L + b, keep_details = FALSE
    )
    for (method in method_names) {
      stability_values[method, b] <- adjusted_rand_index(
        fit_b$labels[[method]], method_consensus[idx, method]
      )
    }
    sub_parts <- do.call(cbind, fit_b$labels)
    sub_grand <- consensus_from_partitions(sub_parts, k = k)$labels
    stability_values["Grand consensus", b] <- adjusted_rand_index(
      sub_grand, grand$labels[idx]
    )
    message(sprintf("  Subsample %d/%d complete.", b, n_subsamples))
  }
  subsample_stability <- rowMeans(stability_values, na.rm = TRUE)

  # Common validation metrics. The true labels are available only for simulation.
  first_distance <- first_details$distance
  all_consensus_labels <- cbind(method_consensus, "Grand consensus" = grand$labels)
  metrics <- do.call(rbind, lapply(colnames(all_consensus_labels), function(method) {
    labels <- all_consensus_labels[, method]
    has_truth <- !all(is.na(truth))
    data.frame(
      method = method,
      ARI = if (has_truth) adjusted_rand_index(labels, truth) else NA_real_,
      NMI = if (has_truth) normalized_mutual_information(labels, truth) else NA_real_,
      best_match_accuracy = if (has_truth) best_label_accuracy(labels, truth) else NA_real_,
      gower_silhouette = mean_silhouette(labels, first_distance),
      imputation_stability = if (method == "Grand consensus") {
        grand_imputation_stability
      } else imputation_stability[method],
      subsample_stability = subsample_stability[method],
      mean_runtime_seconds = if (method == "Grand consensus") {
        sum(rowMeans(runtime_array, na.rm = TRUE))
      } else mean(runtime_array[method, ], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  # Explore K rather than relying on a single criterion. K=4 is used for the
  # benchmark because it is the known simulation truth.
  max_candidate_k <- min(8L, n - 1L)
  message(sprintf("Evaluating candidate cluster counts K=2,...,%d...", max_candidate_k))
  k_grid <- 2:max_candidate_k
  embed <- first_details$embedding$scores
  k_selection <- do.call(rbind, lapply(k_grid, function(k_candidate) {
    set.seed(seed + k_candidate)
    pam_lab <- cluster::pam(first_distance, k = k_candidate, diss = TRUE,
                            cluster.only = TRUE)
    famd_lab <- stats::kmeans(embed, centers = k_candidate,
                              nstart = 20L)$cluster
    gmm_fit <- diag_gmm(embed, k = k_candidate, nstart = 2L,
                        seed = seed + 100L + k_candidate)
    data.frame(
      k = k_candidate,
      gower_pam_silhouette = mean_silhouette(pam_lab, first_distance),
      famd_silhouette = mean(cluster::silhouette(famd_lab, stats::dist(embed))[, "sil_width"]),
      embedding_gmm_bic = gmm_fit$bic
    )
  }))

  # MNAR sensitivity: perturb only originally missing entries on designated
  # features and pool three representative methods over three imputations.
  message("Running directional MNAR sensitivity scenarios...")
  deltas <- c(-1, 0, 1)
  scenario_labels <- matrix(NA_integer_, n, length(deltas),
                            dimnames = list(patient_id, paste0("delta_", deltas)))
  scenario_rows <- list()
  if (sum(dictionary$mnar_direction != 0) == 0L) {
    scenario_labels[,] <- grand$labels
  } else {
    for (s in seq_along(deltas)) {
      delta <- deltas[s]
      scenario_parts <- list()
      for (imp in seq_len(min(3L, m_imputations))) {
        adjusted <- apply_mnar_delta(imputations[[imp]], data, dictionary, delta = delta)
        d <- gower_distance(adjusted, dictionary)
        scenario_parts[[length(scenario_parts) + 1L]] <- cluster::pam(
          d, k = k, diss = TRUE, cluster.only = TRUE
        )
        emb <- mixed_embedding(adjusted, dictionary, max_dim = 15L)$scores
        set.seed(seed + s * 100L + imp)
        scenario_parts[[length(scenario_parts) + 1L]] <- stats::kmeans(
          emb, centers = k, nstart = 20L
        )$cluster
        scenario_parts[[length(scenario_parts) + 1L]] <- latent_class_profile(
          adjusted, dictionary, k = k, nstart = 2L,
          seed = seed + s * 1000L + imp
        )$labels
      }
      lab <- consensus_from_partitions(do.call(cbind, scenario_parts), k = k)$labels
      scenario_labels[, s] <- lab
    }
  }
  baseline_scenario <- scenario_labels[, which(deltas == 0)]
  sensitivity <- do.call(rbind, lapply(seq_along(deltas), function(s) {
    has_truth <- !all(is.na(truth))
    data.frame(
      delta = deltas[s],
      interpretation = c("opposite-direction MNAR", "MAR-style baseline",
                         "designed-direction MNAR")[s],
      ARI_vs_truth = if (has_truth) adjusted_rand_index(scenario_labels[, s], truth) else NA_real_,
      ARI_vs_baseline = adjusted_rand_index(scenario_labels[, s], baseline_scenario),
      reassigned_fraction = 1 - best_label_accuracy(scenario_labels[, s], baseline_scenario),
      stringsAsFactors = FALSE
    )
  }))

  # Descriptive artifacts for the report.
  missingness <- data.frame(
    feature = dictionary$feature,
    type = dictionary$type,
    role = dictionary$role,
    designed_missingness = if ("designed_missingness" %in% names(dictionary)) {
      dictionary$designed_missingness
    } else "unknown",
    missing_rate = vapply(feature_frame(data, dictionary), function(x) mean(is.na(x)), numeric(1L)),
    stringsAsFactors = FALSE
  )
  feature_counts <- as.data.frame(table(dictionary$type, dictionary$role),
                                  stringsAsFactors = FALSE)
  names(feature_counts) <- c("type", "role", "count")
  feature_counts <- feature_counts[feature_counts$count > 0, , drop = FALSE]

  numeric_candidates <- dictionary$feature[
    dictionary$type %in% c("continuous", "discrete", "ordinal")
  ]
  informative_candidates <- dictionary$feature[
    dictionary$type %in% c("continuous", "discrete", "ordinal") &
      dictionary$role == "informative"
  ]
  profile_features <- if (length(informative_candidates)) {
    head(informative_candidates, 18L)
  } else {
    head(numeric_candidates, 18L)
  }
  if (!length(profile_features)) {
    profile_features <- head(dictionary$feature, 18L)
  }
  profile_matrix <- sapply(profile_features, function(nm) {
    x <- imputations[[1L]][[nm]]
    if (is.factor(x)) x <- as.numeric(x)
    as.numeric(scale(as.numeric(x)))
  })
  if (is.null(dim(profile_matrix))) profile_matrix <- matrix(profile_matrix, ncol = 1L)
  profile_matrix[!is.finite(profile_matrix)] <- 0
  profile_means <- t(vapply(sort(unique(grand$labels)), function(g) {
    colMeans(profile_matrix[grand$labels == g, , drop = FALSE])
  }, numeric(ncol(profile_matrix))))
  colnames(profile_means) <- profile_features
  rownames(profile_means) <- paste0("Consensus ", seq_len(nrow(profile_means)))

  certainty <- data.frame(
    patient_id = patient_id,
    consensus_cluster = grand_code,
    membership_certainty = membership_certainty,
    true_cluster = if (!all(is.na(truth))) as.character(truth) else NA_character_,
    stringsAsFactors = FALSE
  )
  confusion <- if (!all(is.na(truth))) table(True = truth, Consensus = grand_code) else NULL

  partition_output <- data.frame(
    patient_id = patient_id,
    true_cluster = if (!all(is.na(truth))) as.character(truth) else NA_character_,
    grand_consensus = grand_code,
    membership_certainty = membership_certainty,
    method_consensus,
    check.names = FALSE
  )

  utils::write.csv(metrics, file.path(results_dir, "method_metrics.csv"), row.names = FALSE)
  utils::write.csv(k_selection, file.path(results_dir, "k_selection.csv"), row.names = FALSE)
  utils::write.csv(sensitivity, file.path(results_dir, "mnar_sensitivity.csv"), row.names = FALSE)
  utils::write.csv(certainty, file.path(results_dir, "membership_certainty.csv"), row.names = FALSE)
  utils::write.csv(partition_output, file.path(results_dir, "cluster_assignments.csv"), row.names = FALSE)
  utils::write.csv(missingness, file.path(results_dir, "missingness_summary.csv"), row.names = FALSE)
  utils::write.csv(as.data.frame(stability_values),
                   file.path(results_dir, "subsample_stability.csv"), row.names = TRUE)

  result <- list(
    settings = list(seed = seed, n = n, p = nrow(dictionary), k = k,
                    m_imputations = m_imputations,
                    n_subsamples = n_subsamples,
                    input_file = normalizePath(data_path, mustWork = FALSE),
                    input_label = basename(data_path),
                    id_column = id_column,
                    truth_column = truth_column,
                    has_truth = !all(is.na(truth)),
                    mnar_feature_count = sum(dictionary$mnar_direction != 0),
                    generated_at = format(Sys.time(), tz = "UTC")),
    data = data,
    first_imputation = imputations[[1L]],
    dictionary = dictionary,
    truth = truth,
    patient_id = patient_id,
    method_names = method_names,
    metrics = metrics,
    k_selection = k_selection,
    sensitivity = sensitivity,
    missingness = missingness,
    feature_counts = feature_counts,
    partitions = partition_array,
    method_consensus = method_consensus,
    method_coassignment = method_coassign,
    grand_consensus = grand_code,
    coassignment = grand$coassignment,
    membership_certainty = membership_certainty,
    stability_values = stability_values,
    profile_means = profile_means,
    confusion = confusion,
    embedding = first_details$embedding,
    runtime_array = runtime_array,
    session_info = utils::capture.output(sessionInfo())
  )
  saveRDS(result, file.path(results_dir, "analysis_results.rds"))
  message("Benchmark complete. Results saved to ", results_dir)
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  data_path <- if (length(args) >= 1L) args[1L] else file.path("data", "mock_patients_typed.rds")
  dictionary_path <- if (length(args) >= 2L) args[2L] else file.path("data", "feature_dictionary.csv")
  run_clustering_benchmark(data_path, dictionary_path)
}
