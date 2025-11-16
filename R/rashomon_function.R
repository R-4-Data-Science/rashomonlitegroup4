# Helper function to compute AIC
.get_aic <- function(model) {
  2 * length(stats::coef(model)) - 2 * as.numeric(stats::logLik(model))
}

# 1. fit_models_dim(X, y, k, m)
fit_models_dim <- function(X, y, k, m) {
  all_vars <- colnames(X)
  if (k > length(all_vars)) {
    return(data.frame())
  }
  
  full_data <- data.frame(y = y, X)
  idx_combs <- utils::combn(seq_along(all_vars), k)
  models_list <- vector("list", ncol(idx_combs))
  
  for (i in seq_len(ncol(idx_combs))) {
    vars_to_use <- all_vars[idx_combs[, i]]
    formula_str <- paste("y ~", paste(vars_to_use, collapse = " + "))
    
    suppressWarnings({
      model <- stats::glm(
        formula = as.formula(formula_str),
        data    = full_data,
        family  = binomial
      )
    })
    
    model_aic     <- .get_aic(model)
    model_summary <- summary(model)
    
    models_list[[i]] <- data.frame(
      vars  = I(list(vars_to_use)),
      aic   = model_aic,
      coef  = I(list(stats::coef(model))),
      pvals = I(list(model_summary$coefficients[, "Pr(>|z|)"]))
    )
  }
  
  models_df <- do.call(rbind, models_list)
  models_df <- models_df[order(models_df$aic), ]
  head(models_df, m)
}

# 2. select_best(models_df, alpha)
select_best <- function(models_df, alpha = 0.5) {
  if (nrow(models_df) == 0) {
    return(data.frame())
  }
  
  n_select <- ceiling(nrow(models_df) * alpha)
  models_df <- models_df[order(models_df$aic), ]
  head(models_df, n_select)
}

# 3. grow_from(best_df, all_vars, k_next, m)
grow_from <- function(best_df, all_vars, k_next, m) {
  if (nrow(best_df) == 0) {
    return(data.frame())
  }
  
  # Use lexical scoping: X and y live in run_rashomon()'s frame
  parent_env <- parent.frame()
  if (!exists("X", envir = parent_env) || !exists("y", envir = parent_env)) {
    stop("grow_from must be called from within run_rashomon(), where X and y are defined.")
  }
  X <- get("X", envir = parent_env)
  y <- get("y", envir = parent_env)
  full_data <- data.frame(y = y, X)
  
  new_models_list <- list()
  
  for (i in seq_len(nrow(best_df))) {
    current_vars   <- unlist(best_df$vars[i])
    available_vars <- setdiff(all_vars, current_vars)
    if (length(available_vars) == 0) next
    
    for (new_var in available_vars) {
      new_vars_set <- sort(c(current_vars, new_var))
      if (length(new_vars_set) != k_next) next  # sanity check
      
      formula_str <- paste("y ~", paste(new_vars_set, collapse = " + "))
      
      suppressWarnings({
        model <- stats::glm(
          formula = as.formula(formula_str),
          data    = full_data,
          family  = binomial
        )
      })
      
      model_aic     <- .get_aic(model)
      model_summary <- summary(model)
      
      new_models_list[[length(new_models_list) + 1L]] <- data.frame(
        vars  = I(list(new_vars_set)),
        aic   = model_aic,
        coef  = I(list(stats::coef(model))),
        pvals = I(list(model_summary$coefficients[, "Pr(>|z|)"]))
      )
    }
  }
  
  if (length(new_models_list) == 0) {
    return(data.frame())
  }
  
  new_models_df <- do.call(rbind, new_models_list)
  
  # Deduplicate models within this dimension
  new_models_df$vars_key <- sapply(
    new_models_df$vars,
    function(v) paste(v, collapse = "|")
  )
  
  best_unique_models_aic <- stats::aggregate(
    aic ~ vars_key,
    data = new_models_df,
    FUN  = min
  )
  
  new_models_df$merge_key <- paste(new_models_df$vars_key, new_models_df$aic, sep = "_")
  best_unique_models_aic$merge_key <- paste(
    best_unique_models_aic$vars_key,
    best_unique_models_aic$aic,
    sep = "_"
  )
  
  unique_indices <- match(best_unique_models_aic$merge_key, new_models_df$merge_key)
  merged_models  <- new_models_df[unique_indices, ]
  
  merged_models$vars_key  <- NULL
  merged_models$merge_key <- NULL
  
  merged_models <- merged_models[order(merged_models$aic), ]
  head(merged_models, m)
}

# Helper functions for run_rashomon summaries

.create_aic_summary <- function(models_df) {
  if (nrow(models_df) == 0) {
    return(data.frame(min = NA, Q1 = NA, median = NA, Q3 = NA, max = NA))
  }
  q <- stats::quantile(models_df$aic, probs = c(0, 0.25, 0.5, 0.75, 1))
  data.frame(min = q[1], Q1 = q[2], median = q[3], Q3 = q[4], max = q[5])
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
  counts
}

# 4. run_rashomon(X, y, pmax, m, alpha)
run_rashomon <- function(X, y, pmax = 5, m = 90, alpha = 0.5) {
  if (pmax < 1) {
    stop("pmax must be at least 1")
  }
  
  all_vars <- colnames(X)
  if (pmax > length(all_vars)) {
    stop("pmax cannot exceed number of predictors in X")
  }
  
  M <- list()
  
  # k = 1: fit all 1-variable models
  k <- 1L
  current_best_df <- fit_models_dim(X, y, k, m)
  selected_df     <- select_best(current_best_df, alpha)
  
  M[[as.character(k)]] <- list(
    best_df          = selected_df,
    model_count      = nrow(selected_df),
    aic_summary      = .create_aic_summary(selected_df),
    predictor_counts = .count_predictors(selected_df, all_vars)
  )
  
  # k = 2..pmax: grow from prior best
  if (pmax >= 2) {
    for (k in 2:pmax) {
      prior_df <- M[[as.character(k - 1L)]]$best_df
      if (nrow(prior_df) == 0) break
      
      current_best_df <- grow_from(
        best_df  = prior_df,
        all_vars = all_vars,
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
  
  M
}
