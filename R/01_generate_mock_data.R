# Generate a mixed-type patient cohort with known clusters and structured missingness.
# Run from the project root with: Rscript R/01_generate_mock_data.R

generate_mock_patients <- function(n = 400L, seed = 20260902L,
                                   output_dir = "data") {
  stopifnot(n >= 100L)
  set.seed(seed)

  k <- 4L
  cluster_sizes <- as.integer(round(n * c(0.29, 0.26, 0.24, 0.21)))
  cluster_sizes[k] <- n - sum(cluster_sizes[-k])
  true_cluster <- factor(rep(seq_len(k), cluster_sizes),
                         levels = seq_len(k),
                         labels = paste0("C", seq_len(k)))
  true_cluster <- sample(true_cluster, n, replace = FALSE)
  cluster_index <- as.integer(true_cluster)

  # Four partially overlapping patient phenotypes in a four-factor space.
  centers <- rbind(
    c( 1.35,  0.15, -0.30,  0.35),
    c( 0.10,  1.45,  0.05, -0.15),
    c(-0.30,  0.40,  1.40,  0.25),
    c(-0.80, -0.75, -0.60, -0.45)
  )
  latent_cor <- matrix(c(
    1.00, 0.25, 0.10, 0.15,
    0.25, 1.00, 0.20, 0.10,
    0.10, 0.20, 1.00, 0.25,
    0.15, 0.10, 0.25, 1.00
  ), 4, 4, byrow = TRUE)
  latent_noise <- matrix(rnorm(n * 4L), n, 4L) %*% chol(latent_cor)
  latent <- centers[cluster_index, , drop = FALSE] + 0.78 * latent_noise
  nuisance <- matrix(rnorm(n * 3L), n, 3L)

  # 30 continuous features: 18 informative, 12 correlated nuisance/noise.
  n_cont <- 30L
  cont <- matrix(NA_real_, n, n_cont)
  for (j in seq_len(n_cont)) {
    if (j <= 18L) {
      loading <- rep(0, 4L)
      loading[(j - 1L) %% 4L + 1L] <- ifelse(j %% 5L == 0L, -0.75, 0.90)
      loading[j %% 4L + 1L] <- 0.25
      score <- as.vector(latent %*% loading) + rnorm(n, sd = 0.95)
    } else {
      score <- 0.55 * nuisance[, (j - 1L) %% 3L + 1L] + rnorm(n, sd = 1.00)
    }
    unit_center <- 25 + 3.5 * j
    unit_scale <- c(1.5, 4, 8, 12, 20)[(j - 1L) %% 5L + 1L]
    cont[, j] <- unit_center + unit_scale * score
  }
  # A few clinically plausible extreme-but-valid measurements.
  outlier_cells <- cbind(sample.int(n, 16L, replace = TRUE),
                         sample.int(n_cont, 16L, replace = TRUE))
  for (r in seq_len(nrow(outlier_cells))) {
    i <- outlier_cells[r, 1L]
    j <- outlier_cells[r, 2L]
    cont[i, j] <- cont[i, j] + sample(c(-1, 1), 1L) *
      4.5 * stats::sd(cont[, j])
  }
  colnames(cont) <- sprintf("continuous_%02d", seq_len(n_cont))

  # 20 overdispersed count features: 10 informative, 10 nuisance/noise.
  n_disc <- 20L
  disc <- matrix(NA_integer_, n, n_disc)
  for (j in seq_len(n_disc)) {
    if (j <= 10L) {
      signal <- 0.42 * latent[, (j - 1L) %% 4L + 1L] +
        0.12 * latent[, j %% 4L + 1L]
    } else {
      signal <- 0.25 * nuisance[, (j - 1L) %% 3L + 1L]
    }
    mu <- pmin(exp(0.45 + 0.065 * j + signal), 30)
    disc[, j] <- stats::rnbinom(n, mu = mu, size = 3.5)
  }
  colnames(disc) <- sprintf("discrete_%02d", seq_len(n_disc))

  # 20 ordered features on a five-level scale: 10 informative, 10 nuisance/noise.
  n_ord <- 20L
  ord <- vector("list", n_ord)
  ordinal_levels <- paste0("L", 1:5)
  for (j in seq_len(n_ord)) {
    if (j <= 10L) {
      score <- 0.82 * latent[, (j - 1L) %% 4L + 1L] + rnorm(n, sd = 1.00)
    } else {
      score <- 0.45 * nuisance[, (j - 1L) %% 3L + 1L] + rnorm(n, sd = 1.10)
    }
    code <- cut(score, breaks = c(-Inf, -1.0, -0.30, 0.30, 1.0, Inf),
                labels = ordinal_levels, ordered_result = TRUE)
    ord[[j]] <- code
  }
  names(ord) <- sprintf("ordinal_%02d", seq_len(n_ord))

  # 20 nominal features with 3-4 categories: 10 informative, 10 noise.
  n_cat <- 20L
  cat_list <- vector("list", n_cat)
  for (j in seq_len(n_cat)) {
    lev <- paste0("K", seq_len(ifelse(j %% 3L == 0L, 3L, 4L)))
    if (j <= 10L) {
      fav <- ((cluster_index + j - 2L) %% length(lev)) + 1L
      p_fav <- 0.56
      sampled <- vapply(seq_len(n), function(i) {
        pr <- rep((1 - p_fav) / (length(lev) - 1L), length(lev))
        pr[fav[i]] <- p_fav
        sample(lev, 1L, prob = pr)
      }, character(1L))
    } else {
      base_prob <- seq_along(lev)
      base_prob <- rev(base_prob / sum(base_prob))
      sampled <- sample(lev, n, replace = TRUE, prob = base_prob)
    }
    cat_list[[j]] <- factor(sampled, levels = lev)
  }
  names(cat_list) <- sprintf("categorical_%02d", seq_len(n_cat))

  # 10 binary features: 5 informative, 5 nuisance/noise.
  n_bin <- 10L
  bin_list <- vector("list", n_bin)
  for (j in seq_len(n_bin)) {
    eta <- if (j <= 5L) {
      -0.35 + 0.95 * latent[, (j - 1L) %% 4L + 1L]
    } else {
      -0.20 + 0.35 * nuisance[, (j - 1L) %% 3L + 1L]
    }
    bin_list[[j]] <- factor(stats::rbinom(n, 1L, stats::plogis(eta)),
                            levels = c(0, 1), labels = c("No", "Yes"))
  }
  names(bin_list) <- sprintf("binary_%02d", seq_len(n_bin))

  complete_features <- data.frame(cont, disc, ord, cat_list, bin_list,
                                  check.names = FALSE)
  feature_names <- names(complete_features)
  feature_type <- c(rep("continuous", n_cont), rep("discrete", n_disc),
                    rep("ordinal", n_ord), rep("categorical", n_cat),
                    rep("binary", n_bin))
  informative <- c(seq_len(n_cont) <= 18L, seq_len(n_disc) <= 10L,
                   seq_len(n_ord) <= 10L, seq_len(n_cat) <= 10L,
                   seq_len(n_bin) <= 5L)
  domain <- sprintf("domain_%02d", ((seq_along(feature_names) - 1L) %% 10L) + 1L)

  # Selected variables have value-dependent observation probabilities. A positive
  # direction means higher unseen values are more likely to be missing.
  mnar_direction <- rep(0L, length(feature_names))
  names(mnar_direction) <- feature_names
  mnar_direction[c("continuous_05", "continuous_14", "continuous_25",
                   "discrete_04", "discrete_13", "ordinal_03",
                   "ordinal_12")] <- c(1L, -1L, 1L, 1L, -1L, 1L, -1L)

  # Site is an observed process variable used only to generate MAR panel effects;
  # it is not included as a clustering feature.
  site <- factor(sample(c("Site_A", "Site_B", "Site_C"), n, replace = TRUE,
                        prob = c(0.42, 0.35, 0.23)))
  care_intensity <- as.numeric(scale(0.55 * latent[, 2L] +
                                     0.35 * latent[, 3L] + rnorm(n)))
  observed_features <- complete_features
  missing_class <- character(length(feature_names))

  for (j in seq_along(feature_names)) {
    x <- complete_features[[j]]
    base_rate <- 0.045 + 0.008 * ((j - 1L) %% 10L)
    eta <- stats::qlogis(base_rate) + 0.40 * care_intensity

    # Panel/site effects create MAR missingness in selected feature families.
    if ((j %% 9L) %in% c(0L, 1L)) {
      eta <- eta + 0.90 * (site == "Site_C") - 0.15 * (site == "Site_A")
    }
    if (mnar_direction[j] != 0L) {
      numeric_x <- if (is.factor(x)) as.numeric(x) else as.numeric(x)
      eta <- eta + 0.90 * mnar_direction[j] * as.numeric(scale(numeric_x))
      missing_class[j] <- "MNAR-sensitive"
    } else if (sd(eta) > 0.05) {
      missing_class[j] <- "MAR-like"
    } else {
      missing_class[j] <- "MCAR-like"
    }
    # A small universal MCAR component, with probabilities capped for usability.
    miss_prob <- pmin(0.60, 0.018 + stats::plogis(eta))
    is_missing <- stats::runif(n) < miss_prob
    observed_features[[j]][is_missing] <- NA
  }

  patient_id <- sprintf("P%04d", seq_len(n))
  observed <- data.frame(patient_id = patient_id,
                         true_cluster = true_cluster,
                         site = site,
                         observed_features,
                         check.names = FALSE)
  complete <- data.frame(patient_id = patient_id,
                         true_cluster = true_cluster,
                         site = site,
                         complete_features,
                         check.names = FALSE)

  dictionary <- data.frame(
    feature = feature_names,
    type = feature_type,
    domain = domain,
    role = ifelse(informative, "informative", "noise_or_nuisance"),
    mnar_direction = unname(mnar_direction),
    designed_missingness = missing_class,
    description = paste("Synthetic", feature_type, "feature"),
    stringsAsFactors = FALSE
  )
  missingness <- data.frame(
    feature = feature_names,
    type = feature_type,
    role = dictionary$role,
    designed_missingness = missing_class,
    missing_rate = vapply(observed_features, function(x) mean(is.na(x)), numeric(1L)),
    stringsAsFactors = FALSE
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(observed, file.path(output_dir, "mock_patients.csv"),
                   row.names = FALSE, na = "")
  utils::write.csv(dictionary, file.path(output_dir, "feature_dictionary.csv"),
                   row.names = FALSE)
  utils::write.csv(missingness, file.path(output_dir, "missingness_summary.csv"),
                   row.names = FALSE)
  saveRDS(observed, file.path(output_dir, "mock_patients_typed.rds"))
  saveRDS(complete, file.path(output_dir, "mock_patients_complete_truth.rds"))

  message(sprintf(
    "Generated %d patients, %d clustering features, 4 clusters, %.1f%% missing cells.",
    n, length(feature_names),
    100 * mean(is.na(observed_features))
  ))
  invisible(list(observed = observed, complete = complete,
                 dictionary = dictionary, missingness = missingness))
}

if (sys.nframe() == 0L) {
  generate_mock_patients()
}
