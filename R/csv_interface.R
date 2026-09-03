# CSV-first interface for the mixed-type clustering benchmark.

source(file.path("R", "02_run_clustering_benchmark.R"))
source(file.path("R", "03_render_report.R"))

csv_help <- function() {
  cat(paste0(
    "Mixed-type clustering from a CSV file\n\n",
    "Usage:\n",
    "  Rscript analyze_csv.R --input PATH [options]\n\n",
    "Required:\n",
    "  --input PATH            CSV with one patient per row\n\n",
    "Common options:\n",
    "  --output-dir PATH       Run directory (default: output/csv_run)\n",
    "  --pdf PATH              PDF path (default: <output-dir>/clustering_report.pdf)\n",
    "  --dictionary PATH       Optional feature dictionary; inferred if omitted\n",
    "  --id-column NAME        Unique patient identifier column\n",
    "  --truth-column NAME     Optional known labels used only for evaluation\n",
    "  --exclude A,B,C         Columns to ignore, such as site or outcomes\n",
    "  --k auto|INTEGER        Number of clusters (default: auto)\n",
    "  --max-k INTEGER         Largest K considered by auto selection (default: 8)\n",
    "  --imputations INTEGER   Completed datasets (default: 8)\n",
    "  --subsamples INTEGER    Repeated 80% stability fits (default: 10)\n",
    "  --seed INTEGER          Reproducible seed (default: 20260902)\n",
    "  --help                  Show this message\n\n",
    "Dictionary columns:\n",
    "  feature, type, domain are required. Allowed types are continuous,\n",
    "  discrete, ordinal, categorical, and binary. Optional columns include\n",
    "  role and mnar_direction.\n"
  ))
}

parse_csv_cli <- function(args) {
  if (!length(args) || any(args %in% c("--help", "-h"))) {
    return(list(help = TRUE))
  }
  known <- c("--input", "--output-dir", "--pdf", "--dictionary",
             "--id-column", "--truth-column", "--exclude", "--k",
             "--max-k", "--imputations", "--subsamples", "--seed")
  out <- list(help = FALSE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (!key %in% known) stop("Unknown option: ", key, ". Use --help.")
    if (i == length(args)) stop("Missing value after ", key)
    value <- args[i + 1L]
    name <- sub("^--", "", key)
    name <- gsub("-", "_", name)
    out[[name]] <- value
    i <- i + 2L
  }
  out
}

infer_id_column <- function(data, requested = NULL) {
  if (!is.null(requested)) {
    if (!requested %in% names(data)) stop("ID column not found: ", requested)
    return(requested)
  }
  common <- c("patient_id", "patientid", "subject_id", "subjectid",
              "record_id", "recordid", "case_id", "caseid", "id")
  normalized <- tolower(gsub("[^a-z0-9]", "", names(data)))
  match_idx <- match(gsub("[^a-z0-9]", "", common), normalized, nomatch = 0L)
  match_idx <- match_idx[match_idx > 0L]
  if (length(match_idx)) return(names(data)[match_idx[1L]])

  high_card_text <- names(data)[vapply(data, function(x) {
    x <- stats::na.omit(x)
    (is.character(x) || is.factor(x)) && length(x) > 0L &&
      length(unique(x)) / length(x) >= 0.98
  }, logical(1L))]
  if (length(high_card_text) == 1L) {
    message("Auto-detected high-cardinality text ID column: ", high_card_text)
    return(high_card_text)
  }
  NULL
}

looks_ordinal <- function(values) {
  values <- tolower(trimws(as.character(values)))
  if (length(values) < 3L || length(values) > 12L) return(FALSE)
  numbered <- grepl("^(l|level|stage|grade|class|score)[ _-]?[0-9]+$", values)
  if (all(numbered)) return(TRUE)
  known_orders <- list(
    c("none", "mild", "moderate", "severe"),
    c("very low", "low", "medium", "high", "very high"),
    c("poor", "fair", "good", "very good", "excellent"),
    c("never", "rarely", "sometimes", "often", "always")
  )
  any(vapply(known_orders, function(ordering) all(values %in% ordering), logical(1L)))
}

infer_feature_type <- function(x) {
  observed <- stats::na.omit(x)
  unique_values <- unique(observed)
  n_unique <- length(unique_values)
  if (is.logical(x) || n_unique <= 2L) return("binary")
  if (is.numeric(x)) {
    integer_like <- all(abs(observed - round(observed)) < 1e-8)
    discrete_limit <- max(30L, ceiling(sqrt(length(observed))))
    if (integer_like && n_unique <= discrete_limit) return("discrete")
    return("continuous")
  }
  if (looks_ordinal(unique_values)) return("ordinal")
  "categorical"
}

prepare_csv_dictionary <- function(data, dictionary_path = NULL,
                                   id_column = NULL, truth_column = NULL,
                                   exclude = character(), output_dir) {
  excluded_metadata <- unique(stats::na.omit(c(id_column, truth_column, exclude)))
  if (!is.null(dictionary_path)) {
    dictionary <- utils::read.csv(dictionary_path, stringsAsFactors = FALSE,
                                  check.names = FALSE)
  } else {
    candidate_names <- setdiff(names(data), excluded_metadata)
    exclusion_rows <- list()
    keep <- logical(length(candidate_names))
    for (i in seq_along(candidate_names)) {
      nm <- candidate_names[i]
      observed <- stats::na.omit(data[[nm]])
      reason <- NULL
      if (!length(observed)) reason <- "all values missing"
      else if (length(unique(observed)) <= 1L) reason <- "constant among observed values"
      if (is.null(reason)) keep[i] <- TRUE
      else exclusion_rows[[length(exclusion_rows) + 1L]] <- data.frame(
        column = nm, reason = reason, stringsAsFactors = FALSE
      )
    }
    features <- candidate_names[keep]
    if (length(features) < 2L) stop("Fewer than two usable feature columns remain.")
    dictionary <- data.frame(
      feature = features,
      type = vapply(data[features], infer_feature_type, character(1L)),
      domain = features,
      role = "unspecified",
      mnar_direction = 0L,
      designed_missingness = "not specified",
      description = "Type inferred from CSV; review before final analysis",
      stringsAsFactors = FALSE
    )
    if (length(exclusion_rows)) {
      utils::write.csv(do.call(rbind, exclusion_rows),
                       file.path(output_dir, "excluded_columns.csv"), row.names = FALSE)
    }
  }

  required <- c("feature", "type", "domain")
  if (!all(required %in% names(dictionary))) {
    stop("Dictionary must contain: ", paste(required, collapse = ", "))
  }
  allowed <- c("continuous", "discrete", "ordinal", "categorical", "binary")
  bad_types <- setdiff(unique(dictionary$type), allowed)
  if (length(bad_types)) stop("Unsupported dictionary type(s): ", paste(bad_types, collapse = ", "))
  missing_features <- setdiff(dictionary$feature, names(data))
  if (length(missing_features)) {
    stop("Dictionary features absent from CSV: ", paste(missing_features, collapse = ", "))
  }
  duplicate_features <- dictionary$feature[duplicated(dictionary$feature)]
  if (length(duplicate_features)) stop("Duplicate dictionary feature: ", duplicate_features[1L])
  if (!"role" %in% names(dictionary)) dictionary$role <- "unspecified"
  if (!"mnar_direction" %in% names(dictionary)) dictionary$mnar_direction <- 0L
  if (!"designed_missingness" %in% names(dictionary)) {
    dictionary$designed_missingness <- "not specified"
  }
  dictionary$mnar_direction[is.na(dictionary$mnar_direction)] <- 0L

  used_path <- file.path(output_dir, "feature_dictionary_used.csv")
  utils::write.csv(dictionary, used_path, row.names = FALSE)
  list(dictionary = dictionary, path = used_path,
       inferred = is.null(dictionary_path))
}

select_k_automatically <- function(data, dictionary, max_k = 8L,
                                   seed = 20260902L) {
  completed <- stochastic_hotdeck_impute(data, dictionary, m = 1L,
                                         seed = seed + 17L)[[1L]]
  d <- gower_distance(completed, dictionary)
  embedding <- mixed_embedding(completed, dictionary, max_dim = 15L)$scores
  candidate_k <- 2:min(max_k, nrow(data) - 1L)
  diagnostics <- do.call(rbind, lapply(candidate_k, function(k) {
    set.seed(seed + k)
    pam_labels <- cluster::pam(d, k = k, diss = TRUE, cluster.only = TRUE)
    latent_labels <- stats::kmeans(embedding, centers = k, nstart = 20L)$cluster
    mixture <- diag_gmm(embedding, k = k, nstart = 2L, seed = seed + 100L + k)
    data.frame(
      k = k,
      gower_pam_silhouette = mean_silhouette(pam_labels, d),
      famd_silhouette = mean(cluster::silhouette(latent_labels,
                                                 stats::dist(embedding))[, "sil_width"]),
      embedding_gmm_bic = mixture$bic
    )
  }))
  diagnostics$aggregate_rank <-
    rank(-diagnostics$gower_pam_silhouette, ties.method = "average") +
    rank(-diagnostics$famd_silhouette, ties.method = "average") +
    rank(diagnostics$embedding_gmm_bic, ties.method = "average")
  selected <- diagnostics$k[which.min(diagnostics$aggregate_rank)]
  list(k = selected, diagnostics = diagnostics)
}

analyze_csv <- function(input_path,
                        output_dir = file.path("output", "csv_run"),
                        pdf_path = NULL,
                        dictionary_path = NULL,
                        id_column = NULL,
                        truth_column = NULL,
                        exclude = character(),
                        k = "auto",
                        max_k = 8L,
                        m_imputations = 8L,
                        n_subsamples = 10L,
                        seed = 20260902L) {
  if (!file.exists(input_path)) stop("Input CSV not found: ", input_path)
  if (!grepl("\\.csv$", input_path, ignore.case = TRUE)) {
    stop("The CSV interface expects a .csv file.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(pdf_path)) pdf_path <- file.path(output_dir, "clustering_report.pdf")
  data <- utils::read.csv(input_path, stringsAsFactors = FALSE,
                          na.strings = c("", "NA", "N/A", "."),
                          check.names = FALSE)
  if (nrow(data) < 20L) stop("At least 20 rows are required for this workflow.")
  id_column <- infer_id_column(data, id_column)
  if (!is.null(truth_column) && !truth_column %in% names(data)) {
    stop("Truth column not found: ", truth_column)
  }
  if (length(exclude)) {
    absent <- setdiff(exclude, names(data))
    if (length(absent)) stop("Excluded column not found: ", absent[1L])
  }

  prepared <- prepare_csv_dictionary(
    data, dictionary_path = dictionary_path, id_column = id_column,
    truth_column = truth_column, exclude = exclude, output_dir = output_dir
  )
  typed <- restore_types(data, prepared$dictionary)

  if (identical(tolower(as.character(k)), "auto")) {
    message("Selecting K automatically from three diagnostics...")
    auto <- select_k_automatically(typed, prepared$dictionary,
                                   max_k = as.integer(max_k), seed = seed)
    selected_k <- auto$k
    utils::write.csv(auto$diagnostics,
                     file.path(output_dir, "auto_k_diagnostics.csv"), row.names = FALSE)
    message("Automatic K selected: ", selected_k)
  } else {
    selected_k <- as.integer(k)
    if (!is.finite(selected_k) || selected_k < 2L || selected_k >= nrow(data)) {
      stop("--k must be 'auto' or an integer between 2 and nrow(data)-1.")
    }
    stale_auto_path <- file.path(output_dir, "auto_k_diagnostics.csv")
    if (file.exists(stale_auto_path)) file.remove(stale_auto_path)
  }

  results_dir <- file.path(output_dir, "results")
  result <- run_clustering_benchmark(
    data_path = input_path,
    dictionary_path = prepared$path,
    results_dir = results_dir,
    seed = as.integer(seed),
    k = selected_k,
    m_imputations = as.integer(m_imputations),
    n_subsamples = as.integer(n_subsamples),
    id_column = id_column,
    truth_column = truth_column
  )
  result$settings$dictionary_inferred <- prepared$inferred
  result$settings$dictionary_path <- normalizePath(prepared$path, mustWork = FALSE)
  saveRDS(result, file.path(results_dir, "analysis_results.rds"))
  render_benchmark_report(file.path(results_dir, "analysis_results.rds"), pdf_path)

  manifest <- data.frame(
    item = c("input_csv", "feature_dictionary", "results_directory", "pdf_report",
             "selected_k", "id_column", "truth_column"),
    value = c(normalizePath(input_path, mustWork = FALSE),
              normalizePath(prepared$path, mustWork = FALSE),
              normalizePath(results_dir, mustWork = FALSE),
              normalizePath(pdf_path, mustWork = FALSE),
              selected_k,
              ifelse(is.null(id_column), "", id_column),
              ifelse(is.null(truth_column), "", truth_column)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(output_dir, "run_manifest.csv"), row.names = FALSE)
  message("CSV analysis complete. PDF: ", pdf_path)
  invisible(result)
}

run_csv_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- parse_csv_cli(args)
  if (isTRUE(opts$help)) {
    csv_help()
    return(invisible(NULL))
  }
  if (is.null(opts$input)) stop("--input is required. Use --help for usage.")
  exclude <- if (!is.null(opts$exclude) && nzchar(opts$exclude)) {
    trimws(strsplit(opts$exclude, ",", fixed = TRUE)[[1L]])
  } else character()
  analyze_csv(
    input_path = opts$input,
    output_dir = opts$output_dir %||% file.path("output", "csv_run"),
    pdf_path = opts$pdf %||% NULL,
    dictionary_path = opts$dictionary %||% NULL,
    id_column = opts$id_column %||% NULL,
    truth_column = opts$truth_column %||% NULL,
    exclude = exclude,
    k = opts$k %||% "auto",
    max_k = as.integer(opts$max_k %||% 8L),
    m_imputations = as.integer(opts$imputations %||% 8L),
    n_subsamples = as.integer(opts$subsamples %||% 10L),
    seed = as.integer(opts$seed %||% 20260902L)
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
