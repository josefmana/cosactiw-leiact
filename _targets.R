library(targets)
library(tarchetypes)

# Load packages
tar_option_set(packages = c(
  "here",        # path listing
  "tidyverse",   # data wrangling
  "openxlsx",    # data reading
  "effsize",     # effect sizes
  "rcompanion",  # computation of some effect sizes
  "performance", # model diagnostics
  "brms",        # Bayesian regressions
  "bayestestR",  # ETIs
  "bayesplot",   # plotting
  "patchwork"    # arranging plot
))

# Load functions
tar_source()

# Use all cores for model fitting
options(mc.cores = parallel::detectCores())

# Build the pipeline
list(
  target_data,
  target_counts,
  target_intensities,
  tar_quarto(report, "report.qmd", quiet = FALSE)
)
