# Helper function to compute AIC
.get_aic <- function(model) {
  return(2 * length(stats::coef(model)) - 2 * as.numeric(stats::logLik(model)))
}

# ---
# 1. fit_models_dim(x, y, k, m)
#' @export
fit_models_dim <- function(x, y, k, m) {
  all_vars <- colnames(x)
  if (k > length(all_vars)) {
    return(data.frame())
  }

  # Get all possible combinations of k predictors
  idx_combs <- utils::combn(seq_along(all_vars), k)

  models_list <- list()

  # Loop through all combinations and fit the models
  for (i in 1:ncol(idx_combs)) {
    vars_to_use <- all_vars[idx_combs[, i]]

    # Create the formula string for glm
    formula_str <- paste("y ~", paste(vars_to_use, collapse = " + "))

    suppressWarnings({
      model <- stats::glm(as.formula(formula_str),
                          data = data.frame(x),
                          family = binomial)
    })

    model_aic <- .get_aic(model)
    model_summary <- summary(model)

    # Store results, using I() to keep lists intact in the data frame column
    models_list[[i]] <- data.frame(
      vars = I(list(vars_to_use)),
      aic = model_aic,
      coef = I(list(stats::coef(model))),
      pvals = I(list(model_summary$coefficients[, "Pr(>|z|)"]))
    )
  }

  models_df <- do.call(rbind, models_list)

  # Order by AIC (lower is better) and return the top m
  models_df <- models_df[order(models_df$aic), ]

  return(head(models_df, m))
}

# 2. select_best(models_df, alpha)
#' @export
select_best <- function(models_df, alpha = 0.5) {
  if (nrow(models_df) == 0) {
    return(data.frame())
  }

  # Calculate number of models to select (top alpha% by count)
  n_select <- ceiling(nrow(models_df) * alpha)

  # Order by AIC (lower is better) and return top alpha%
  models_df <- models_df[order(models_df$aic), ]

  return(head(models_df, n_select))
}

# 3. grow_from(best_df, all_vars, k_next, m)
#' @export
grow_from <- function(best_df, all_vars, k_next, m) {
  if (nrow(best_df) == 0) {
    return(data.frame())
  }

  new_models_list <- list()

  # I. Generate and Fit New Models by "growing" the prior best set
  for (i in 1:nrow(best_df)) {
    current_vars <- unlist(best_df$vars[i])
    available_vars <- setdiff(all_vars, current_vars)

    for (new_var in available_vars) {
      new_vars_set <- sort(c(current_vars, new_var)) # Sorted for uniqueness check
      formula_str <- paste("y ~", paste(new_vars_set, collapse = " + "))

      suppressWarnings({
        model <- stats::glm(as.formula(formula_str),
                            data = data.frame(x),
                            family = binomial)
      })

      model_aic <- .get_aic(model)
      model_summary <- summary(model)

      new_models_list[[length(new_models_list) + 1]] <- data.frame(
        vars = I(list(new_vars_set)),
        aic = model_aic,
        coef = I(list(stats::coef(model))),
        pvals = I(list(model_summary$coefficients[, "Pr(>|z|)"]))
      )
    }
  }

  new_models_df <- do.call(rbind, new_models_list)

  if (nrow(new_models_df) == 0) {
    return(data.frame())
  }

  # II. Prevent Duplicates (Crucial Fix)
  new_models_df$vars_key <- sapply(new_models_df$vars,
                                   function(v) paste(v, collapse = "|"))

  # Find the minimum AIC for each unique combination of variables
  best_unique_models_aic <- stats::aggregate(aic ~ vars_key,
                                             data = new_models_df,
                                             FUN = min)

  # Create a composite key for robust filtering
  new_models_df$merge_key <- paste(new_models_df$vars_key, new_models_df$aic, sep="_")
  best_unique_models_aic$merge_key <- paste(best_unique_models_aic$vars_key, best_unique_models_aic$aic, sep="_")

  # Use match() to find the index of the single, best instance of each unique model
  unique_indices <- match(best_unique_models_aic$merge_key, new_models_df$merge_key)
  merged_models <- new_models_df[unique_indices, ]

  # III. Clean up and select top m
  merged_models$vars_key <- NULL
  merged_models$merge_key <- NULL

  merged_models <- merged_models[order(merged_models$aic), ]

  return(head(merged_models, m))
}

# Helper functions for run_rashomon summaries

.create_aic_summary <- function(models_df) {
  if (nrow(models_df) == 0) {
    return(data.frame(min = NA, Q1 = NA, median = NA, Q3 = NA, max = NA))
  }
  q <- stats::quantile(models_df$aic, probs = c(0, 0.25, 0.5, 0.75, 1))
  return(data.frame(min = q[1], Q1 = q[2], median = q[3], Q3 = q[4], max = q[5]))
}

.count_predictors <- function(models_df, all_vars) {
  counts <- stats::setNames(rep(0, length(all_vars)), all_vars)

  if (nrow(models_df) == 0) {
    return(counts)
  }

  for (vars_list in models_df$vars) {
    for (v in unlist(vars_list)) {
      if (v %in% names(counts)) {
        counts[v] <- counts[v] + 1
      }
    }
  }
  return(counts)
}

# 4. run_rashomon(x, y, pmax, s, m, alpha)
run_rashomon <- function(x, y, pmax = 5, s = 1, m = 90, alpha = 0.5) {
  all_vars <- colnames(x)
  M <- list()
  k <- s
  
  # dimension k = s
  current_best_df <- fit_models_dim(x, y, k, m)
  selected_df     <- select_best(current_best_df, alpha)
  
  M[[as.character(k)]] <- list(
    best_df          = selected_df,
    model_count      = nrow(selected_df),
    aic_summary      = .create_aic_summary(selected_df),
    predictor_counts = .count_predictors(selected_df, all_vars)
  )
  
  # grow to dimensions s+1, ..., pmax
  if (pmax > s) {
    for (k in (s + 1):pmax) {
      prior_df <- M[[as.character(k - 1)]]$best_df
      
      current_best_df <- grow_from(
        best_df  = prior_df,
        all_vars = all_vars,
        x        = x,
        y        = y,
        k_next   = k,
        m        = m
      )
      
      selected_df <- select_best(current_best_df, alpha)
      
      M[[as.character(k)]] <- list(
        best_df          = selected_df,
        model_count      = nrow(selected_df),
        aic_summary      = .create_aic_summary(selected_df),
        predictor_counts = .count_predictors(selected_df, all_vars)
      )
      
      if (nrow(selected_df) == 0) break
    }
  }
  return(M)
}
