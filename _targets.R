library(targets)
library(tarchetypes)

# Load packages
tar_option_set(packages = c(
  "here",        # path listing
  "tidyverse",   # data wrangling
  "rlang",       # tidy evaluation
  "openxlsx",    # data reading
  "apa",         # APA-formatted output
  "effsize",     # effect sizes
  "rcompanion",  # computation of some effect sizes
  "performance", # model diagnostics
  "brms",        # Bayesian regressions
  "bayestestR",  # ETIs
  "bayesplot",   # plotting
  "patchwork",   # arranging plot
  "gt",          # APA-formatted tables
  "gtExtras"     # gt extensions
))

# Load functions
tar_source()

# Use all cores for model fitting
options(mc.cores = parallel::detectCores())

# Build the pipeline
list(
  target_files,
  target_data,
  target_counts,
  target_intensities,
  target_cobra,
  tar_quarto(report, "outputs/report.qmd", quiet = FALSE)
)
