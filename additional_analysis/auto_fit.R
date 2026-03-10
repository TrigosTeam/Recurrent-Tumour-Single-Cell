library(glmmTMB)
library(dplyr)

# ============================================================
# TRUE BULLETPROOF safe_fit — intercepts C++ level errors
# ============================================================

safe_fit <- function(formula_str, data) {
  
  fit_attempt <- function(formula_str, data, family_use, optimizer_name) {
    
    # Redirect ALL output including C++ stderr
    fit <- NULL
    
    invisible(capture.output({
      fit <- try(
        suppressWarnings(
          suppressMessages(
            switch(optimizer_name,
                   "BFGS" = glmmTMB(
                     as.formula(formula_str), data = data,
                     family = family_use, zi = ~0, dispformula = ~1,
                     control = glmmTMBControl(
                       optimizer = optim, optArgs = list(method = "BFGS"))
                   ),
                   "NM" = glmmTMB(
                     as.formula(formula_str), data = data,
                     family = family_use, zi = ~0, dispformula = ~1,
                     control = glmmTMBControl(
                       optimizer = optim, optArgs = list(method = "Nelder-Mead"))
                   ),
                   "nlminb" = glmmTMB(
                     as.formula(formula_str), data = data,
                     family = family_use, zi = ~0, dispformula = ~1,
                     control = glmmTMBControl(optimizer = nlminb)
                   )
            )
          )
        ),
        silent = TRUE   # ← this is key: silences C++ level errors
      )
    }, type = "message"))  # ← capture stderr too
    
    # Return NULL if try() caught an error
    if (inherits(fit, "try-error")) return(NULL)
    
    # Validate: AIC must be finite
    aic_val <- try(AIC(fit), silent = TRUE)
    if (inherits(aic_val, "try-error")) return(NULL)
    if (is.na(aic_val) || !is.finite(aic_val)) return(NULL)
    
    return(fit)
  }
  
  # ── Attempt matrix: optimizer x family ───────────────────
  combos <- list(
    list(family = betabinomial(), opt = "BFGS"),
    list(family = betabinomial(), opt = "NM"),
    list(family = betabinomial(), opt = "nlminb"),
    list(family = binomial(),     opt = "BFGS"),
    list(family = binomial(),     opt = "NM"),
    list(family = binomial(),     opt = "nlminb")
  )
  
  for (combo in combos) {
    fit <- fit_attempt(formula_str, data,
                       family_use    = combo$family,
                       optimizer_name = combo$opt)
    if (!is.null(fit)) return(fit)
  }
  
  return(NULL)  # all 6 attempts failed
}


# ============================================================
# auto_select_model (robust version)
# ============================================================

auto_select_model <- function(data, response, fixed_var = "pathology") {
  
  candidates <- list(
    site_only = paste0(response, " ~ site"),
    fixed_only          = paste0(response, " ~ ", fixed_var),
    fixed_site          = paste0(response, " ~ ", fixed_var, " + site"),
    fixed_patient       = paste0(response, " ~ ", fixed_var, " + patient"),
    fixed_both          = paste0(response, " ~ ", fixed_var, " + patient + site"),
    re_site             = paste0(response, " ~ ", fixed_var, " + (1 | site)"),
    re_patient          = paste0(response, " ~ ", fixed_var, " + (1 | patient)"),
    re_patient_site     = paste0(response, " ~ ", fixed_var, " + (1 | patient/site)"),
    re_site_fix_patient = paste0(response, " ~ ", fixed_var, " + patient + (1 | site)"),
    re_patient_fix_site = paste0(response, " ~ ", fixed_var, " + site + (1 | patient)")
  )
  
  cat("Fitting", length(candidates), "candidate models...\n")
  
  fits       <- list()
  fit_status <- list()
  
  for (nm in names(candidates)) {
    cat(" ", nm, "... ")
    fit <- safe_fit(candidates[[nm]], data)
    if (!is.null(fit)) {
      fits[[nm]]       <- fit
      fit_status[[nm]] <- "OK"
      cat("✓\n")
    } else {
      fit_status[[nm]] <- "FAILED"
      cat("✗ (all optimizers failed)\n")
    }
  }
  
  # ── Require at least 1 successful fit ─────────────────────
  if (length(fits) == 0) {
    cat("\n⛔ No models converged for this cell type. Skipping.\n")
    return(NULL)
  }
  
  # ── Metrics table ─────────────────────────────────────────
  metrics_df <- lapply(names(fits), function(nm) {
    m <- fits[[nm]]
    data.frame(
      model   = nm,
      formula = deparse(formula(m)),
      df      = attr(logLik(m), "df"),
      AIC     = AIC(m),
      BIC     = BIC(m),
      logLik  = as.numeric(logLik(m)),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    arrange(AIC) %>%
    mutate(
      delta_AIC  = AIC - min(AIC),
      delta_BIC  = BIC - min(BIC),
      AIC_weight = round(exp(-0.5 * delta_AIC) /
                           sum(exp(-0.5 * delta_AIC)) * 100, 1)
    )
  
  best_name  <- metrics_df$model[1]
  best_model <- fits[[best_name]]
  
  # ── Print results ─────────────────────────────────────────
  cat("\n📊 MODEL COMPARISON\n")
  print(
    metrics_df %>%
      select(model, df, AIC, BIC, delta_AIC, AIC_weight) %>%
      mutate(across(where(is.numeric), ~round(., 3))),
    row.names = FALSE
  )
  
  cat("\n✅ BEST MODEL:", best_name)
  cat("\n   Formula   :", deparse(formula(best_model)))
  cat("\n   AIC       :", round(AIC(best_model), 3))
  cat("\n   BIC       :", round(BIC(best_model), 3), "\n")
  
  # ── Decision rationale ────────────────────────────────────
  cat("\n💡 DECISION RATIONALE\n")
  if (nrow(metrics_df) > 1 && !is.na(metrics_df$delta_AIC[2])) {
    d <- metrics_df$delta_AIC[2]
    if      (d > 4) cat("  → Strong support for best model (ΔAIC =",   round(d,2), ")\n")
    else if (d > 2) cat("  → Moderate support for best model (ΔAIC =", round(d,2), ")\n")
    else            cat("  → Models competitive (ΔAIC =", round(d,2),
                        ") — preferring simpler model\n")
  } else {
    cat("  → Only one model converged\n")
  }
  
  invisible(list(
    best_model = best_model,
    best_name  = best_name,
    metrics    = metrics_df,
    all_fits   = fits,
    status     = fit_status
  ))
}


# ============================================================
# MAIN LOOP
# ============================================================

anno_list <- unique(df$cell_anno)
results   <- list()

for (anno in anno_list) {
  cat("\n", strrep("═", 55), "\n")
  cat("  CELL TYPE:", anno, "\n")
  cat(strrep("═", 55), "\n")
  
  df_sub <- df %>% filter(cell_anno == anno)
  
  if (nrow(df_sub) < 6 || n_distinct(df_sub$pathology) < 2) {
    cat("⚠ Skipping — insufficient data\n")
    next
  }
  
  res <- auto_select_model(
    data     = df_sub,
    response = "cbind(n, total_normal - n)"
  )
  
  if (!is.null(res)) results[[anno]] <- res
}

# ── Final summary table ───────────────────────────────────
summary_table <- lapply(names(results), function(anno) {
  m <- results[[anno]]$metrics[1, ]
  data.frame(
    cell_anno  = anno,
    best_model = results[[anno]]$best_name,
    AIC        = round(m$AIC, 2),
    BIC        = round(m$BIC, 2),
    AIC_weight = paste0(m$AIC_weight, "%")
  )
}) %>% bind_rows()

cat("\n\n📌 FINAL SUMMARY\n")
cat(strrep("=", 60), "\n")
print(summary_table, row.names = FALSE)
