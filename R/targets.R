# --- Read and format data -----------------------------------------------------
target_data <- list(
  targets::tar_target(
    name    = mapping_file,
    command = here("data-raw", "COSACTIW_DCERA-aktivity2.xlsx"),
    format  = "file"
  ),
  targets::tar_target(
    name    = data_file,
    command = here("data-raw", "COSACTIW_NANOK_pro-jamovi.xlsx"),
    format  = "file"
  ),
  targets::tar_target(
    name    = sa_file,
    command = here("data-raw", "info.rds"),
    format  = "file"
  ),
  targets::tar_target(
    name    = data_long,
    command = get_data(
      mfile = mapping_file,
      dfile = data_file,
      sfile = sa_file
    )
  ),
  targets::tar_target(
    name    = data_wide,
    command = pivot_data(data_long)
  ),
  targets::tar_target(
    name    = preprocessed_data,
    command = preprocess_data(data_long)
  ),
  targets::tar_target(
    name    = ppc_data,
    command = add_ppc_categories(preprocessed_data)
  )
)

# --- Analyse activity counts --------------------------------------------------
target_counts <- list(
  targets::tar_target(
    name    = N,
    command = count_subjects(data_wide)
  ),
  targets::tar_target(
    name    = demography_comparisons,
    command = compare_demography(data_wide)
  ),
  targets::tar_target(
    name    = activity_counts,
    command = count_activities(
      input = data_long,
      data  = data_wide
    )
  ),
  targets::tar_target(
    name    = counts_comparisons_table,
    command = add_comparisons(
      table = activity_counts,
      data  = data_wide
    )
  ),
  targets::tar_target(
    name    = counts_regression_specifications,
    command = specify_regression(data_long)
  ),
  targets::tar_target(
    name    = counts_regressions,
    command = count_regressions(
      data  = data_wide,
      specs = counts_regression_specifications
    )
  ),
  targets::tar_target(
    name    = counts_regressions_diagnostics,
    command = diagnose_regressions(counts_regressions)
  ),
  targets::tar_target(
    name    = counts_regressions_table,
    command = add_regression_parameters(
      table = activity_counts,
      fits  = counts_regressions,
      input = data_long
    )
  ),
  targets::tar_target(
    name    = cross_tables,
    command = contingency_tables(
      input = data_long,
      data  = data_wide
    )
  ),
  targets::tar_target(
    name    = chisquares,
    command = chi_squares(cross_tables)
  )
)

# --- Analyse activity intensities ---------------------------------------------
target_intensities <- list(
  targets::tar_target(
    name    = priors,
    command = set_priors()
  ),
  targets::tar_target(
    name    = formulas,
    command = set_formulas()
  ),
  targets::tar_target(
    name    = intensities_loglinear_regressions,
    command = fit_regressions(
      formulas = formulas,
      priors   = priors,
      data     = preprocessed_data
    )
  ),
  targets::tar_target(
    name    = loo_comparisons,
    command = psis_loo(intensities_loglinear_regressions)
  ),
  targets::tar_target(
    name    = posterior_predictive_checks,
    command = perform_posterior_checks(
      fits = intensities_loglinear_regressions,
      data = ppc_data,
      save = TRUE
    )
  ),
  targets::tar_target(
    name    = posterior_expectations,
    command = compute_posterior_expectations(
      data   = preprocessed_data,
      fit    = intensities_loglinear_regressions$ordered_time,
      output = "expectations"
    )
  ),
  targets::tar_target(
    name    = posterior_contrasts,
    command = compute_posterior_expectations(
      data   = preprocessed_data,
      fit    = intensities_loglinear_regressions$ordered_time,
      output = "contrasts"
    )
  ),
  targets::tar_target(
    name    = posterior_interaction_plots,
    command = draw_interaction_plots(
      posterior_expect = posterior_expectations,
      text = FALSE,
      save = TRUE
    )
  )
)
