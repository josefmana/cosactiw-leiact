library(targets)
library(tarchetypes)

tar_option_set(packages = c(
  "here",        # for path listing
  "tidyverse",   # for data wrangling
  "openxlsx",    # for data reading
  "effsize",     # for effect sizes
  "rcompanion",  # for computation of some effect sizes
  "performance", # for model diagnostics
  "brms",        # for Bayesian regressions
  "bayestestR",  # for ETIs
  "bayesplot",   # for plotting
  "patchwork"    # for arranging plot
))

tar_source()

options(mc.cores = parallel::detectCores())

list(
  # READ & FORMAT DATA ----
  tar_target(
    name    = mapping_file,
    command = here("data-raw", "COSACTIW_DCERA-aktivity2.xlsx"),
    format  = "file"
  ),
  tar_target(
    name    = data_file,
    command = here("data-raw", "COSACTIW_NANOK_pro-jamovi.xlsx"),
    format  = "file"
  ),
  tar_target(
    name    = sa_file,
    command = here("data-raw", "info.rds"),
    format  = "file"
  ),
  tar_target(
    name    = data_long,
    command = get_data(
      mfile = mapping_file,
      dfile = data_file,
      sfile = sa_file
    )
  ),
  tar_target(
    name    = data_wide,
    command = pivot_data(data_long)
  ),

  # ACTIVITY COUNTS ----
  tar_target(
    name    = N,
    command = count_subjects(data_wide)
  ),
  tar_target(
    name    = demography_comparisons,
    command = compare_demography(data_wide)
  ),
  tar_target(
    name    = activity_counts,
    command = count_activities(
      input = data_long,
      data  = data_wide
    )
  ),
  tar_target(
    name    = counts_comparisons_table,
    command = add_comparisons(
      table = activity_counts,
      data  = data_wide
    )
  ),
  tar_target(
    name    = counts_regression_specifications,
    command = specify_regression(data_long)
  ),
  tar_target(
    name    = counts_regressions,
    command = count_regressions(
      data  = data_wide,
      specs = counts_regression_specifications
    )
  ),
  tar_target(
    name    = counts_regressions_diagnostics,
    command = diagnose_regressions(counts_regressions)
  ),
  tar_target(
    name    = counts_regressions_table,
    command = add_regression_parameters(
      table = activity_counts,
      fits  = counts_regressions,
      input = data_long
    )
  ),
  tar_target(
    name    = cross_tables,
    command = contingency_tables(
      input = data_long,
      data  = data_wide
    )
  ),
  tar_target(
    name    = chisquares,
    command = chi_squares(cross_tables)
  ),

  # INTENSITIES ----
  tar_target(
    name    = preprocessed_data,
    command = preprocess_data(data_long)
  ),
  tar_target(
    name    = priors,
    command = set_priors()
  ),
  tar_target(
    name    = formulas,
    command = set_formulas()
  ),
  tar_target(
    name    = intensities_loglinear_regressions,
    command = fit_regressions(
      formulas = formulas,
      priors   = priors,
      data     = preprocessed_data
    )
  ),
  tar_target(
    name    = loo_comparisons,
    command = psis_loo(intensities_loglinear_regressions)
  ),
  tar_target(
    name    = ppc_data,
    command = add_ppc_categories(preprocessed_data)
  ),
  tar_target(
    name    = posterior_predictive_checks,
    command = perform_posterior_checks(
      fits = intensities_loglinear_regressions,
      data = ppc_data,
      save = TRUE
    )
  ),
  tar_target(
    name    = posterior_expectations,
    command = compute_posterior_expectations(
      data   = preprocessed_data,
      fit    = intensities_loglinear_regressions$ordered_time,
      output = "expectations"
    )
  ),
  tar_target(
    name    = posterior_contrasts,
    command = compute_posterior_expectations(
      data   = preprocessed_data,
      fit    = intensities_loglinear_regressions$ordered_time,
      output = "contrasts"
    )
  ),
  tar_target(
    name    = posterior_interaction_plots,
    command = draw_interaction_plots(
      posterior_expect = posterior_expectations,
      text = FALSE,
      save = TRUE
    )
  ),

  # REPORT ----
  tar_quarto(
    path  = "report.qmd",
    name  = report,
    quiet = FALSE
  )
)
