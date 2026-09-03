# Dependency-light functions for clustering mixed-type data.
# Only the recommended R package `cluster` is required.

require_cluster_package <- function() {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("The recommended R package 'cluster' is required.")
  }
}

restore_types <- function(data, dictionary) {
  for (i in seq_len(nrow(dictionary))) {
    nm <- dictionary$feature[i]
    tp <- dictionary$type[i]
    if (!nm %in% names(data)) stop("Missing feature in data: ", nm)
    if (tp %in% c("continuous", "discrete")) {
      data[[nm]] <- as.numeric(data[[nm]])
    } else if (tp == "ordinal") {
      vals <- unique(stats::na.omit(as.character(data[[nm]])))
      suppressWarnings(num_vals <- as.numeric(sub("^[^0-9.-]*", "", vals)))
      lev <- if (all(is.finite(num_vals))) vals[order(num_vals)] else sort(vals)
      data[[nm]] <- ordered(as.character(data[[nm]]), levels = lev)
    } else {
      data[[nm]] <- factor(data[[nm]])
    }
  }
  data
}

feature_frame <- function(data, dictionary) {
  data[, dictionary$feature, drop = FALSE]
}

domain_weights <- function(dictionary) {
  counts <- table(dictionary$domain)
  w <- 1 / as.numeric(counts[dictionary$domain])
  w / mean(w)
}

gower_distance <- function(data, dictionary) {
  require_cluster_package()
  x <- feature_frame(data, dictionary)
  x <- restore_types(x, dictionary)
  suppressWarnings(cluster::daisy(
    x,
    metric = "gower",
    weights = domain_weights(dictionary)
  ))
}

stochastic_hotdeck_impute <- function(data, dictionary, m = 8L,
                                      donors = 10L, seed = 1L) {
  stopifnot(m >= 1L, donors >= 2L)
  set.seed(seed)
  x <- restore_types(data, dictionary)
  d <- as.matrix(gower_distance(x, dictionary))
  diag(d) <- Inf
  output <- vector("list", m)

  for (imp in seq_len(m)) {
    completed <- x
    for (j in seq_len(nrow(dictionary))) {
      nm <- dictionary$feature[j]
      miss <- which(is.na(x[[nm]]))
      if (!length(miss)) next
      obs <- which(!is.na(x[[nm]]))
      if (!length(obs)) stop("Feature is completely missing: ", nm)

      donor_index <- vapply(miss, function(i) {
        ord <- obs[order(d[i, obs], na.last = NA)]
        if (!length(ord)) ord <- obs
        cand <- head(ord, min(donors, length(ord)))
        local_d <- d[i, cand]
        local_d[!is.finite(local_d)] <- max(local_d[is.finite(local_d)], 1)
        temperature <- max(stats::median(local_d[is.finite(local_d)]), 0.03)
        pr <- exp(-local_d / temperature)
        sample(cand, 1L, prob = pr)
      }, integer(1L))

      completed[[nm]][miss] <- x[[nm]][donor_index]
      if (dictionary$type[j] == "continuous") {
        spread <- stats::sd(x[[nm]], na.rm = TRUE)
        if (is.finite(spread) && spread > 0) {
          completed[[nm]][miss] <- completed[[nm]][miss] +
            stats::rnorm(length(miss), sd = 0.025 * spread)
        }
      }
    }
    if (anyNA(feature_frame(completed, dictionary))) {
      stop("Internal error: imputation left missing values.")
    }
    output[[imp]] <- completed
  }
  output
}

numeric_block <- function(data, dictionary, include_ordinal = TRUE) {
  types <- c("continuous", "discrete", if (include_ordinal) "ordinal")
  dict <- dictionary[dictionary$type %in% types, , drop = FALSE]
  mats <- lapply(seq_len(nrow(dict)), function(i) {
    x <- data[[dict$feature[i]]]
    if (dict$type[i] == "discrete") x <- log1p(as.numeric(x))
    else if (dict$type[i] == "ordinal") x <- as.numeric(x)
    else x <- as.numeric(x)
    s <- stats::sd(x)
    if (!is.finite(s) || s < 1e-8) rep(0, length(x)) else (x - mean(x)) / s
  })
  out <- do.call(cbind, mats)
  colnames(out) <- dict$feature
  out
}

mixed_embedding <- function(data, dictionary, max_dim = 15L,
                            variance_target = 0.75) {
  # A dependency-light FAMD/PCAmix-style embedding. Numeric and ordinal variables
  # are standardized. Each nominal variable is represented by centered/scaled
  # indicators and normalized so variables with many levels do not dominate.
  num <- numeric_block(data, dictionary, include_ordinal = TRUE)
  cat_dict <- dictionary[dictionary$type %in% c("categorical", "binary"), , drop = FALSE]
  dummy_blocks <- vector("list", nrow(cat_dict))
  if (nrow(cat_dict)) {
    for (j in seq_len(nrow(cat_dict))) {
      f <- droplevels(factor(data[[cat_dict$feature[j]]]))
      mm <- stats::model.matrix(~ f - 1)
      p <- colMeans(mm)
      mm <- sweep(mm, 2L, p, "-")
      mm <- sweep(mm, 2L, sqrt(pmax(p * (1 - p), 0.02)), "/")
      mm <- mm / sqrt(max(1, ncol(mm) - 1L))
      colnames(mm) <- paste0(cat_dict$feature[j], "::", levels(f))
      dummy_blocks[[j]] <- mm
    }
  }
  design <- cbind(num, if (length(dummy_blocks)) do.call(cbind, dummy_blocks))
  design[!is.finite(design)] <- 0
  fit <- stats::prcomp(design, center = FALSE, scale. = FALSE,
                       rank. = min(max_dim, nrow(design) - 1L, ncol(design)))
  variance <- fit$sdev^2 / sum(fit$sdev^2)
  target_dim <- which(cumsum(variance) >= variance_target)[1L]
  target_dim <- min(max_dim, max(5L, target_dim), ncol(fit$x))
  list(scores = fit$x[, seq_len(target_dim), drop = FALSE],
       full_scores = fit$x,
       variance = variance,
       dimensions = target_dim,
       loadings = fit$rotation,
       design = design)
}

mode_value <- function(x) {
  tab <- table(x, useNA = "no")
  names(tab)[which.max(tab)]
}

kprototypes_cluster <- function(data, dictionary, k, nstart = 12L,
                                max_iter = 80L, lambda = 1, seed = 1L) {
  set.seed(seed)
  num <- numeric_block(data, dictionary, include_ordinal = TRUE)
  cat_names <- dictionary$feature[dictionary$type %in% c("categorical", "binary")]
  cats <- lapply(cat_names, function(nm) as.character(data[[nm]]))
  names(cats) <- cat_names
  n <- nrow(data)
  best <- NULL

  for (start in seq_len(nstart)) {
    init <- sample.int(n, k)
    num_centers <- num[init, , drop = FALSE]
    cat_centers <- if (length(cats)) {
      lapply(cats, function(x) x[init])
    } else list()
    labels <- rep(NA_integer_, n)

    for (iter in seq_len(max_iter)) {
      cost_num <- sapply(seq_len(k), function(g) {
        rowMeans((sweep(num, 2L, num_centers[g, ], "-"))^2)
      })
      if (!is.matrix(cost_num)) cost_num <- matrix(cost_num, ncol = k)
      if (length(cats)) {
        cost_cat <- sapply(seq_len(k), function(g) {
          mismatch <- vapply(seq_len(n), function(i) {
            mean(vapply(seq_along(cats), function(j) {
              cats[[j]][i] != cat_centers[[j]][g]
            }, logical(1L)))
          }, numeric(1L))
        })
      } else {
        cost_cat <- matrix(0, n, k)
      }
      new_labels <- max.col(-(cost_num + lambda * cost_cat), ties.method = "random")
      if (identical(new_labels, labels)) break
      labels <- new_labels

      for (g in seq_len(k)) {
        members <- which(labels == g)
        if (!length(members)) {
          members <- sample.int(n, 1L)
          labels[members] <- g
        }
        num_centers[g, ] <- colMeans(num[members, , drop = FALSE])
        if (length(cats)) {
          for (j in seq_along(cats)) {
            cat_centers[[j]][g] <- mode_value(cats[[j]][members])
          }
        }
      }
    }
    total_cost <- sum(cost_num[cbind(seq_len(n), labels)] +
                        lambda * cost_cat[cbind(seq_len(n), labels)])
    if (is.null(best) || total_cost < best$objective) {
      best <- list(labels = labels, objective = total_cost,
                   num_centers = num_centers, cat_centers = cat_centers)
    }
  }
  best
}

logsumexp_rows <- function(x) {
  mx <- apply(x, 1L, max)
  mx + log(rowSums(exp(x - mx)))
}

diag_gmm <- function(x, k, nstart = 4L, max_iter = 120L,
                     tol = 1e-5, seed = 1L) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  best <- NULL
  set.seed(seed)
  global_var <- pmax(apply(x, 2L, stats::var), 0.15)

  for (start in seq_len(nstart)) {
    km <- stats::kmeans(x, centers = k, nstart = 1L,
                        iter.max = 50L)
    z <- matrix(0, n, k)
    z[cbind(seq_len(n), km$cluster)] <- 1
    old_ll <- -Inf

    for (iter in seq_len(max_iter)) {
      nk <- pmax(colSums(z), 1e-8)
      pi_g <- nk / n
      mu <- t(z) %*% x / nk
      variances <- matrix(NA_real_, k, p)
      for (g in seq_len(k)) {
        centered <- sweep(x, 2L, mu[g, ], "-")
        variances[g, ] <- colSums(centered^2 * z[, g]) / nk[g]
      }
      variances <- pmax(variances, matrix(0.08 * global_var, k, p, byrow = TRUE))

      logp <- matrix(NA_real_, n, k)
      for (g in seq_len(k)) {
        centered <- sweep(x, 2L, mu[g, ], "-")
        logp[, g] <- log(pmax(pi_g[g], 1e-10)) -
          0.5 * rowSums(log(2 * pi * variances[g, ]) +
                          centered^2 / variances[g, ])
      }
      norm <- logsumexp_rows(logp)
      ll <- sum(norm)
      z <- exp(logp - norm)
      if (is.finite(old_ll) && abs(ll - old_ll) < tol * (1 + abs(old_ll))) break
      old_ll <- ll
    }
    npar <- (k - 1L) + 2L * k * p
    fit <- list(labels = max.col(z), posterior = z, loglik = ll,
                bic = -2 * ll + npar * log(n), means = mu,
                variances = variances, iterations = iter)
    if (is.null(best) || fit$loglik > best$loglik) best <- fit
  }
  best
}

latent_class_profile <- function(data, dictionary, k, nstart = 3L,
                                 max_iter = 100L, tol = 1e-5,
                                 alpha = 0.5, seed = 1L) {
  # Conditional-independence finite mixture: Gaussian profiles for standardized
  # continuous/log-count features and multinomials for ordinal/nominal/binary data.
  num <- numeric_block(data, dictionary, include_ordinal = FALSE)
  cat_names <- dictionary$feature[dictionary$type %in%
                                    c("ordinal", "categorical", "binary")]
  cat_codes <- lapply(cat_names, function(nm) as.integer(factor(data[[nm]])))
  nlevels_cat <- vapply(cat_names, function(nm) nlevels(factor(data[[nm]])), integer(1L))
  n <- nrow(data)
  p <- ncol(num)
  global_var <- pmax(apply(num, 2L, stats::var), 0.15)
  best <- NULL
  set.seed(seed)
  init_embed <- mixed_embedding(data, dictionary, max_dim = 10L)$scores

  for (start in seq_len(nstart)) {
    initial <- stats::kmeans(init_embed, centers = k, nstart = 1L)$cluster
    z <- matrix(0, n, k)
    z[cbind(seq_len(n), initial)] <- 1
    old_ll <- -Inf

    for (iter in seq_len(max_iter)) {
      nk <- pmax(colSums(z), 1e-8)
      pi_g <- nk / n
      mu <- t(z) %*% num / nk
      variances <- matrix(NA_real_, k, p)
      for (g in seq_len(k)) {
        centered <- sweep(num, 2L, mu[g, ], "-")
        variances[g, ] <- colSums(centered^2 * z[, g]) / nk[g]
      }
      variances <- pmax(variances, matrix(0.08 * global_var, k, p, byrow = TRUE))
      theta <- vector("list", length(cat_codes))
      for (j in seq_along(cat_codes)) {
        onehot <- stats::model.matrix(~ factor(cat_codes[[j]],
                                               levels = seq_len(nlevels_cat[j])) - 1)
        probs <- t(z) %*% onehot + alpha
        probs <- probs / rowSums(probs)
        theta[[j]] <- probs
      }

      logp <- matrix(NA_real_, n, k)
      for (g in seq_len(k)) {
        centered <- sweep(num, 2L, mu[g, ], "-")
        lp <- log(pmax(pi_g[g], 1e-10)) -
          0.5 * rowSums(log(2 * pi * variances[g, ]) +
                          centered^2 / variances[g, ])
        for (j in seq_along(cat_codes)) {
          lp <- lp + log(pmax(theta[[j]][g, cat_codes[[j]]], 1e-12))
        }
        logp[, g] <- lp
      }
      norm <- logsumexp_rows(logp)
      ll <- sum(norm)
      z <- exp(logp - norm)
      if (is.finite(old_ll) && abs(ll - old_ll) < tol * (1 + abs(old_ll))) break
      old_ll <- ll
    }
    category_parameters <- sum(nlevels_cat - 1L)
    npar <- (k - 1L) + 2L * k * p + k * category_parameters
    fit <- list(labels = max.col(z), posterior = z, loglik = ll,
                bic = -2 * ll + npar * log(n), iterations = iter)
    if (is.null(best) || fit$loglik > best$loglik) best <- fit
  }
  best
}

spectral_gower <- function(dissimilarity, k, seed = 1L) {
  d <- as.matrix(dissimilarity)
  off_diag <- d[row(d) != col(d)]
  sigma <- stats::median(off_diag[is.finite(off_diag) & off_diag > 0])
  sigma <- max(sigma, 1e-4)
  affinity <- exp(-(d^2) / (2 * sigma^2))
  diag(affinity) <- 0
  degree <- pmax(rowSums(affinity), 1e-12)
  normalized <- affinity / sqrt(outer(degree, degree))
  eig <- eigen(normalized, symmetric = TRUE)
  u <- eig$vectors[, seq_len(k), drop = FALSE]
  u <- u / pmax(sqrt(rowSums(u^2)), 1e-12)
  set.seed(seed)
  labels <- stats::kmeans(u, centers = k, nstart = 25L, iter.max = 100L)$cluster
  list(labels = labels, embedding = u, sigma = sigma,
       eigenvalues = eig$values[seq_len(k)])
}

run_all_methods <- function(data, dictionary, k, seed = 1L,
                            keep_details = FALSE) {
  method_names <- c("Gower + PAM", "Gower + hierarchical", "k-prototypes",
                    "FAMD-style + k-means", "Latent-embedding mixture",
                    "Latent class/profile", "Gower spectral")
  labels <- setNames(vector("list", length(method_names)), method_names)
  runtimes <- setNames(numeric(length(method_names)), method_names)
  details <- list()

  timing <- function(name, expression) {
    start <- proc.time()[[3L]]
    value <- force(expression)
    runtimes[name] <<- proc.time()[[3L]] - start
    value
  }

  d <- timing("Gower + PAM", gower_distance(data, dictionary))
  labels[["Gower + PAM"]] <- timing(
    "Gower + PAM",
    cluster::pam(d, k = k, diss = TRUE, cluster.only = TRUE)
  )
  # Include the shared distance computation in both distance-method runtimes.
  gower_time <- attr(d, "timing")
  hc <- timing("Gower + hierarchical", stats::hclust(d, method = "average"))
  labels[["Gower + hierarchical"]] <- stats::cutree(hc, k = k)

  kp <- timing("k-prototypes", kprototypes_cluster(
    data, dictionary, k = k, seed = seed + 11L
  ))
  labels[["k-prototypes"]] <- kp$labels

  embed <- timing("FAMD-style + k-means", mixed_embedding(data, dictionary))
  set.seed(seed + 21L)
  km <- timing("FAMD-style + k-means", stats::kmeans(
    embed$scores, centers = k, nstart = 30L, iter.max = 100L
  ))
  labels[["FAMD-style + k-means"]] <- km$cluster

  gmm <- timing("Latent-embedding mixture", diag_gmm(
    embed$scores, k = k, seed = seed + 31L
  ))
  labels[["Latent-embedding mixture"]] <- gmm$labels

  lcp <- timing("Latent class/profile", latent_class_profile(
    data, dictionary, k = k, seed = seed + 41L
  ))
  labels[["Latent class/profile"]] <- lcp$labels

  spec <- timing("Gower spectral", spectral_gower(d, k = k, seed = seed + 51L))
  labels[["Gower spectral"]] <- spec$labels

  if (keep_details) {
    details <- list(distance = d, hierarchical = hc, kprototypes = kp,
                    embedding = embed, gmm = gmm, latent_class_profile = lcp,
                    spectral = spec)
  }
  list(labels = labels, runtimes = runtimes, details = details)
}

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2
  n <- sum(tab)
  index <- sum(choose2(tab))
  expected <- sum(choose2(rowSums(tab))) * sum(choose2(colSums(tab))) / choose2(n)
  maximum <- 0.5 * (sum(choose2(rowSums(tab))) + sum(choose2(colSums(tab))))
  if (maximum == expected) return(1)
  (index - expected) / (maximum - expected)
}

normalized_mutual_information <- function(a, b) {
  tab <- table(a, b)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(vapply(seq_len(nrow(nz)), function(i) {
    r <- nz[i, 1L]; c <- nz[i, 2L]
    pxy[r, c] * log(pxy[r, c] / (px[r] * py[c]))
  }, numeric(1L)))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(0)
  mi / sqrt(hx * hy)
}

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    rest <- all_permutations(x[-i])
    cbind(x[i], rest)
  }))
}

align_labels <- function(predicted, reference) {
  pred_text <- as.character(predicted)
  ref_text <- as.character(reference)
  pred_levels <- sort(unique(pred_text))
  ref_levels <- sort(unique(ref_text))
  if (length(pred_levels) != length(ref_levels) || length(pred_levels) > 8L) {
    return(as.integer(factor(predicted)))
  }
  pred_code <- as.integer(factor(pred_text, levels = pred_levels))
  ref_code <- as.integer(factor(ref_text, levels = ref_levels))
  perms <- all_permutations(seq_along(ref_levels))
  scores <- apply(perms, 1L, function(p) sum(p[pred_code] == ref_code))
  best <- perms[which.max(scores), ]
  best[pred_code]
}

best_label_accuracy <- function(predicted, truth) {
  tab <- table(factor(predicted), factor(truth))
  dimension <- max(dim(tab))
  padded <- matrix(0, dimension, dimension)
  padded[seq_len(nrow(tab)), seq_len(ncol(tab))] <- tab
  permutations <- all_permutations(seq_len(dimension))
  matched <- apply(permutations, 1L, function(permutation) {
    sum(padded[cbind(seq_len(dimension), permutation)])
  })
  max(matched) / sum(tab)
}

consensus_from_partitions <- function(partitions, k) {
  partitions <- as.matrix(partitions)
  n <- nrow(partitions)
  coassignment <- matrix(0, n, n)
  for (j in seq_len(ncol(partitions))) {
    lab <- partitions[, j]
    coassignment <- coassignment + outer(lab, lab, "==")
  }
  coassignment <- coassignment / ncol(partitions)
  diag(coassignment) <- 1
  labels <- cluster::pam(stats::as.dist(1 - coassignment), k = k,
                         diss = TRUE, cluster.only = TRUE)
  list(labels = labels, coassignment = coassignment)
}

pairwise_partition_ari <- function(partitions) {
  partitions <- as.matrix(partitions)
  if (ncol(partitions) < 2L) return(NA_real_)
  pairs <- utils::combn(seq_len(ncol(partitions)), 2L)
  mean(apply(pairs, 2L, function(idx) {
    adjusted_rand_index(partitions[, idx[1L]], partitions[, idx[2L]])
  }))
}

mean_silhouette <- function(labels, dissimilarity) {
  if (length(unique(labels)) < 2L) return(NA_real_)
  mean(cluster::silhouette(as.integer(factor(labels)), dissimilarity)[, "sil_width"])
}

select_k_automatically <- function(data, dictionary, max_k = 8L,
                                   seed = 20260902L) {
  max_k <- as.integer(max_k)
  if (!is.finite(max_k) || max_k < 2L) stop("max_k must be at least 2.")
  completed <- stochastic_hotdeck_impute(data, dictionary, m = 1L,
                                         seed = seed + 17L)[[1L]]
  d <- gower_distance(completed, dictionary)
  embedding <- mixed_embedding(completed, dictionary, max_dim = 15L)$scores
  candidate_k <- 2:min(max_k, nrow(data) - 1L)
  diagnostics <- do.call(rbind, lapply(candidate_k, function(k) {
    set.seed(seed + k)
    pam_labels <- cluster::pam(d, k = k, diss = TRUE, cluster.only = TRUE)
    latent_labels <- stats::kmeans(embedding, centers = k, nstart = 20L)$cluster
    mixture <- diag_gmm(embedding, k = k, nstart = 2L,
                        seed = seed + 100L + k)
    data.frame(
      k = k,
      gower_pam_silhouette = mean_silhouette(pam_labels, d),
      famd_silhouette = mean(cluster::silhouette(
        latent_labels, stats::dist(embedding))[, "sil_width"]),
      embedding_gmm_bic = mixture$bic
    )
  }))
  diagnostics$aggregate_rank <-
    rank(-diagnostics$gower_pam_silhouette, ties.method = "average") +
    rank(-diagnostics$famd_silhouette, ties.method = "average") +
    rank(diagnostics$embedding_gmm_bic, ties.method = "average")
  best_rank <- min(diagnostics$aggregate_rank)
  selected <- min(diagnostics$k[diagnostics$aggregate_rank == best_rank])
  list(k = selected, diagnostics = diagnostics)
}

apply_mnar_delta <- function(imputed, original, dictionary, delta = 0.75) {
  out <- imputed
  targets <- which(dictionary$mnar_direction != 0)
  for (j in targets) {
    nm <- dictionary$feature[j]
    miss <- which(is.na(original[[nm]]))
    if (!length(miss)) next
    direction <- dictionary$mnar_direction[j]
    tp <- dictionary$type[j]
    if (tp == "continuous") {
      spread <- stats::sd(original[[nm]], na.rm = TRUE)
      out[[nm]][miss] <- out[[nm]][miss] + delta * direction * spread
    } else if (tp == "discrete") {
      spread <- stats::sd(original[[nm]], na.rm = TRUE)
      out[[nm]][miss] <- pmax(0, round(out[[nm]][miss] +
                                        delta * direction * spread))
    } else if (tp == "ordinal") {
      lev <- levels(out[[nm]])
      code <- as.integer(out[[nm]][miss]) + round(delta * direction)
      code <- pmax(1L, pmin(length(lev), code))
      out[[nm]][miss] <- ordered(lev[code], levels = lev)
    }
  }
  out
}
