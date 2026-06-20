#' Pre-process activity-intensity data for Bayesian regression
#'
#' @description
#' Takes the long-format activity records from the raw data list, joins
#' demographic and cognitive covariates, and applies the variable
#' transformations required before fitting the Bayesian mixed-effects models
#' in [fit_regressions()]. Key operations include:
#' \itemize{
#'   \item Recoding `Intensity` to a 1–5 integer scale (dropping `"do not
#'     know"`) and computing `logIntensity`.
#'   \item Deriving `Time_past` (years before assessment) alongside ordered
#'     and numeric representations of `Time_bin`.
#'   \item Creating composite labels for posterior predictive check grouping
#'     (`Type_Season_Time`, `SA_Type_Time`, `SA_Season_Time`).
#'   \item Setting sum-to-zero contrasts (±0.5) for `SA`, `Seasonal`, and
#'     `Activity_type`.
#' }
#'
#' @param input A named list as returned by [get_data()], containing at
#'   minimum the elements `act` (long-format activity data), `dem`
#'   (demographic data with `Age_years`), and `cog` (cognitive data with `SA`
#'   and other covariates).
#'
#' @return A data frame in long format (one row per participant × activity ×
#'   time-bin observation) with the following additional columns beyond the
#'   raw `act` data:
#'   \describe{
#'     \item{`Age_years`}{Numeric. Age in years (joined from `dem`).}
#'     \item{`Intensity`}{Integer 1–5 (or `NA` for `"do not know"`), shifted
#'       so the scale begins at 1.}
#'     \item{`logIntensity`}{Numeric. Natural log of `Intensity`.}
#'     \item{`Time_past`}{Numeric. `Time_bin - Age_years` (years before
#'       assessment).}
#'     \item{`Time_bin`}{Ordered factor with levels from `act$Time_bin`.}
#'     \item{`Time_num`}{Integer. Numeric index of `Time_bin`.}
#'     \item{`Type_Season_Time`, `SA_Type_Time`, `SA_Season_Time`}{Character.
#'       Composite grouping labels for PPCs.}
#'   }
#'   Sum-to-zero contrasts (±0.5) are set on `SA`, `Seasonal`, and
#'   `Activity_type`.
#'
#' @importFrom dplyr left_join filter mutate
#'
#' @seealso [set_priors()], [set_formulas()], [fit_regressions()]
#'
#' @export
preprocess_data <- function(input) {

  with(input, act %>%

         left_join(dem[, c("ID", "Age_years")]) %>%
         filter(complete.cases(Intensity)) %>%
         left_join(cog) %>%
         mutate(

           # re-format categorical variables
           Seasonal      = factor(Seasonal),
           Activity_type = factor(Activity_type),
           Activity      = factor(Activity),
           ID            = factor(ID),

           # intensity/frequency (outcome) — shift so scale starts at 1
           Intensity    = as.integer(if_else(Intensity == "do not know", NA, Intensity)) - 1,
           logIntensity = log(as.numeric(Intensity)),

           # time (predictor) variables
           Time_past = Time_bin - Age_years,
           Time_bin  = factor(Time_bin, levels = unique(act$Time_bin), ordered = TRUE),
           Time_num  = as.numeric(Time_bin),

           # PPC grouping labels
           Type_Season_Time = paste(Activity_type, Seasonal, Time_bin, sep = "_"),
           SA_Type_Time     = paste(SA, Activity_type, Time_bin, sep = "_"),
           SA_Season_Time   = paste(SA, Seasonal, Time_bin, sep = "_")

         ) %>%

         # set sum-to-zero contrasts (nonSA = -0.5, SA = 0.5; non-seasonal = -0.5, etc.)
         within(., {
           contrasts(SA)            <- -contr.sum(2) / 2
           contrasts(Seasonal)      <- -contr.sum(2) / 2
           contrasts(Activity_type) <- -contr.sum(2) / 2
         }))
}


#' Construct weakly informative priors for the intensity models
#'
#' @description
#' Returns a \pkg{brms} prior specification intended for the Bayesian
#' mixed-effects models fit by [fit_regressions()]. Priors are weakly
#' informative on the log scale of intensity:
#' \describe{
#'   \item{Intercept}{Normal(1, 3)}
#'   \item{Fixed effects (`b`)}{Normal(0, 3)}
#'   \item{Random-effect SDs (`sd`)}{Exponential(2)}
#'   \item{Residual SD (`sigma`)}{Exponential(2)}
#'   \item{Correlation matrices (`cor`)}{LKJ(2)}
#' }
#'
#' @return A `brmsprior` data frame (as returned by [brms::prior()]) that can
#'   be passed directly to the `prior` argument of [brms::brm()].
#'
#' @importFrom brms prior
#'
#' @seealso [fit_regressions()], [set_formulas()]
#'
#' @export
set_priors <- function() {
  c(
    prior(normal(1, 3),   class = "Intercept"),
    prior(normal(0, 3),   class = "b"),
    prior(exponential(2), class = "sd"),
    prior(exponential(2), class = "sigma"),
    prior(lkj(2),         class = "cor")
  )
}


#' Define model formulas for the intensity mixed-effects models
#'
#' @description
#' Returns a named list of \pkg{brms} formula objects representing three
#' parameterisations of time for the log-intensity model:
#' \describe{
#'   \item{`linear_time`}{Time bin entered as a numeric linear predictor
#'     (`Time_num`).}
#'   \item{`discrete_time`}{Time entered as years before assessment
#'     (`Time_past`), so that the slope reflects change relative to the
#'     participant's age.}
#'   \item{`ordered_time`}{Time bin entered as an ordered monotonic predictor
#'     (`mo(Time_bin)`) using \pkg{brms}'s isotonic-regression prior.}
#' }
#' All formulas include a three-way interaction with SuperAging status (`SA`)
#' and activity type (`Activity_type`), plus by-participant and by-category
#' random slopes for the time predictor.
#'
#' @return A named list of `brmsformula` objects.
#'
#' @importFrom brms bf mo
#'
#' @seealso [fit_regressions()], [set_priors()]
#'
#' @export
set_formulas <- function() {
  list(
    linear_time   = bf(logIntensity ~ 1 +    Time_num  * SA * Activity_type + (1 +    Time_num  | ID) + (1 +    Time_num  | Category)),
    discrete_time = bf(logIntensity ~ 1 +    Time_bin  * SA * Activity_type + (1 +    Time_num  | ID) + (1 +    Time_num  | Category)),
    ordered_time  = bf(logIntensity ~ 1 + mo(Time_bin) * SA * Activity_type + (1 + mo(Time_bin) | ID) + (1 + mo(Time_bin) | Category))
  )
}


#' Fit a set of Bayesian mixed-effects models via \pkg{brms}
#'
#' @description
#' Iterates over a named list of \pkg{brms} formula objects (as returned by
#' [set_formulas()]) and fits each as a Gaussian mixed-effects model with the
#' supplied prior. A fixed random seed (`87542`) is used for reproducibility.
#'
#' @param data A pre-processed data frame as returned by [preprocess_data()].
#' @param priors A `brmsprior` object as returned by [set_priors()].
#' @param formulas A named list of `brmsformula` objects as returned by
#'   [set_formulas()].
#'
#' @return A named list of `brmsfit` objects, one per formula, sharing the
#'   same names as `formulas`.
#'
#' @importFrom brms brm gaussian
#'
#' @seealso [set_formulas()], [set_priors()], [psis_loo()]
#'
#' @export
fit_regressions <- function(data, priors, formulas) {
  lapply(
    set_names(names(formulas)),
    function(i) brm(
      formula = formulas[[i]],
      family  = gaussian(link = "identity"),
      prior   = priors,
      data    = data,
      seed    = 87542,
      file    = paste0("models/", i)
    )
  )
}


#' Compare intensity models via PSIS-LOO cross-validation
#'
#' @description
#' Applies \pkg{loo}'s [loo::loo()] to the three fitted models produced by
#' [fit_regressions()] (`linear_time`, `discrete_time`, and `ordered_time`) to
#' compute Pareto-smoothed importance-sampling leave-one-out (PSIS-LOO)
#' estimates and rank the models by their expected log predictive density
#' (ELPD).
#'
#' @param fits A named list of `brmsfit` objects as returned by
#'   [fit_regressions()], expected to contain at minimum the elements
#'   `linear_time`, `discrete_time`, and `ordered_time`.
#'
#' @return A `loolist` comparison object as returned by [loo::loo()].
#'
#' @importFrom loo loo
#'
#' @seealso [fit_regressions()]
#'
#' @export
psis_loo <- function(fits) {
  with(fits, loo(linear_time, discrete_time, ordered_time))
}


#' Add posterior predictive check grouping labels to a data frame
#'
#' @description
#' Appends three composite character columns to a data frame for use as
#' grouping variables in posterior predictive checks (PPCs):
#' \describe{
#'   \item{`Type_Season_Time`}{Combination of activity type, seasonality, and
#'     time bin.}
#'   \item{`SA_Type_Time`}{Combination of SuperAging status, activity type,
#'     and time bin.}
#'   \item{`SA_Season_Time`}{Combination of SuperAging status, seasonality,
#'     and time bin.}
#' }
#'
#' @param data A data frame containing the columns `Activity_type`,
#'   `Seasonal`, `Time_bin`, and `SA`.
#'
#' @return `data` with the three new character columns appended.
#'
#' @seealso [perform_posterior_checks()]
#'
#' @export
add_ppc_categories <- function(data) {
  mutate(
    data,
    Type_Season_Time = paste(Activity_type, Seasonal, Time_bin, sep = "_"),
    SA_Type_Time     = paste(SA, Activity_type, Time_bin, sep = "_"),
    SA_Season_Time   = paste(SA, Seasonal, Time_bin, sep = "_")
  )
}


#' Perform and optionally save posterior predictive check plots
#'
#' @description
#' For each fitted model in `fits`, produces three PPC plots:
#' \describe{
#'   \item{`dens`}{Density-overlay plot comparing observed and predicted
#'     distributions of `logIntensity`, grouped by `SA_Type_Time`.}
#'   \item{`mean`}{Grouped statistic plot for the mean.}
#'   \item{`sd`}{Grouped statistic plot for the standard deviation.}
#' }
#' Plots are produced via [ppc_plot()] and, when `save = TRUE`, written as
#' 300 dpi JPEG files to `"figures/"` (created automatically if absent).
#'
#' @param fits A named list of `brmsfit` objects as returned by
#'   [fit_regressions()].
#' @param data A data frame as returned by [preprocess_data()] (or with PPC
#'   categories added via [add_ppc_categories()]), containing `logIntensity`
#'   and `SA_Type_Time`.
#' @param save Logical. If `TRUE` (default), each plot is saved to
#'   `"figures/ppc_<stat>_<model>.jpg"`.
#'
#' @return A nested named list: outer names are model names (from `fits`),
#'   inner names are `"dens"`, `"mean"`, and `"sd"`. Each leaf is a `ggplot`
#'   object.
#'
#' @importFrom ggplot2 ggsave last_plot
#' @importFrom here here
#'
#' @seealso [ppc_plot()], [fit_regressions()], [add_ppc_categories()]
#'
#' @export
perform_posterior_checks <- function(fits, data, save = TRUE) {

  lapply(
    set_names(names(fits)),
    function(y) lapply(
      set_names(x = c("dens", "mean", "sd")),
      function(f) {

        # build the plot
        if (f == "dens") {
          plt <- ppc_plot(
            fit   = fits[[y]],
            data  = subset(data, complete.cases(logIntensity)),
            y     = "logIntensity",
            x     = "SA_Type_Time",
            labs  = list(
              "Observed (thick lines) and predicted (thin lines) distributions of log(Intensity) of recalled leisure activities",
              "Cells represent different combinations of SA status, activity type, and time bin."
            ),
            meth  = "dens_overlay_grouped",
            draws = sample(1:4e3, 1e2),
            stat  = NULL
          )
        } else {
          plt <- ppc_plot(
            fit   = fits[[y]],
            data  = subset(data, complete.cases(logIntensity)),
            y     = "logIntensity",
            x     = "SA_Type_Time",
            labs  = list(
              paste0("Observed (thick bars) and predicted (histograms) ", f, " of log(Intensity) of recalled leisure activities"),
              "Cells represent different combinations of SA status, activity type, and time bin."
            ),
            meth  = "stat_grouped",
            draws = 1:4e3,
            stat  = f
          )
        }

        # optionally save
        if (save) {
          new_folder("figures")
          ggsave(
            filename = here("figures", paste0("ppc_", f, "_", y, ".jpg")),
            plot     = plt,
            dpi      = 300,
            width    = 12.6,
            height   = 13.3
          )
        }

        return(plt)
      }
    )
  )
}


#' Compute posterior expectations and pairwise contrasts
#'
#' @description
#' Extracts population-level posterior predictions from a fitted model for a
#' grid of `Activity_type`, `SA`, and `Time_bin` combinations (with random
#' effects marginalised out). Then summarises either:
#' \describe{
#'   \item{`expectations`}{Median and 95% ETI for each cell on the log and
#'     back-transformed (exp) scales.}
#'   \item{`contrasts`}{Pairwise differences (and their two-way interactions)
#'     between levels of `Activity_type`, `SA`, and/or `Time_bin`, on both
#'     scales. Each contrast is summarised as median [ETI] and the probability
#'     of direction (pd).}
#' }
#'
#' @param data A data frame as returned by [preprocess_data()], used to
#'   determine factor levels for the prediction grid.
#' @param fit A single `brmsfit` object as returned by [fit_regressions()].
#' @param output A character string: either `"expectations"` (default) or
#'   `"contrasts"`. Determines which object is returned.
#'
#' @return
#' When `output = "expectations"`: a data frame in long format with columns
#' `Activity_type`, `SA`, `Time_bin`, `scale` (`"log"` or `"explog"`),
#' `Estimate`, `ETI_low`, and `ETI_high`.
#'
#' When `output = "contrasts"`: a data frame with one row per contrast × scale
#' combination and columns `diff`, `diff2`, `timebin`, `actype`, `SA`,
#' `scale`, `Estimate_log`, `pd_log`, `Estimate_logexp`, `pd_logexp`.
#'
#' @importFrom brms posterior_epred
#' @importFrom bayestestR eti p_direction
#' @importFrom tidyr pivot_longer pivot_wider expand_grid
#' @importFrom dplyr mutate
#'
#' @seealso [fit_regressions()], [draw_interaction_plots()]
#'
#' @export
compute_posterior_expectations <- function(data, fit, output = "expectations") {

  ## ---- Prediction grid ----
  d_seq <- with(
    data,
    expand.grid(
      ID            = NA,
      Category      = NA,
      Activity_type = c(levels(Activity_type), NA),
      SA            = c(levels(SA), NA),
      Time_bin      = levels(Time_bin)
    )
  ) %>%
    mutate(Time_bin = as.ordered(Time_bin)) %>%
    within(., {
      contrasts(SA)            <- -contr.sum(2) / 2
      contrasts(Activity_type) <- -contr.sum(2) / 2
    })

  # draw population-level posterior predictions
  ppred <- posterior_epred(
    object     = fit,
    newdata    = d_seq,
    re_formula = NA
  ) %>%
    `colnames<-`(with(d_seq, paste(SA, Activity_type, Time_bin, sep = "_")))

  ## ---- Posterior expectations ----
  expectations <- cbind(
    d_seq,
    apply(
      X      = ppred,
      MARGIN = 2,
      FUN    = function(x) c(
        logxEstimate    = median(x),
        logxETI_low     = quantile(x, prob = .025, names = FALSE),
        logxETI_high    = quantile(x, prob = .975, names = FALSE),
        explogxEstimate = median(exp(x)),
        explogxETI_low  = quantile(exp(x), prob = .025, names = FALSE),
        explogxETI_high = quantile(exp(x), prob = .975, names = FALSE)
      )
    ) %>% t()
  ) %>%
    as.data.frame() %>%
    pivot_longer(
      cols          = starts_with("log") | starts_with("explog"),
      names_to      = c("scale", "term"),
      names_pattern = "(.*)x(.*)",
      values_to     = "value"
    ) %>%
    pivot_wider(
      names_from  = "term",
      values_from = "value",
      id_cols     = c("Activity_type", "SA", "Time_bin", "scale")
    )

  ## ---- Pairwise contrasts specification ----
  post_comp <- list(

    # main effects and simple effects of type and SA
    expand_grid(
      diff    = c("physical - mental", "SA - nonSA", "physical_SA - mental_SA", "physical_nonSA - mental_nonSA"),
      diff2   = NA,
      timebin = levels(d_seq$Time_bin),
      actype  = NA,
      SA      = NA
    ),

    # type × SA interaction
    expand_grid(
      diff    = "physical - mental",
      diff2   = "SA - nonSA",
      timebin = levels(d_seq$Time_bin),
      actype  = NA,
      SA      = NA
    ),

    # time contrasts
    expand_grid(
      diff = data.frame(v2 = seq(30, 80, 5)) %>%
        mutate(v1 = v2 + 5, diff = paste0(v1, " - ", v2)) %>%
        select(diff) %>%
        unlist(use.names = FALSE) %>%
        c("85 - 30"),
      diff2   = NA,
      timebin = NA,
      actype  = c("mental", "physical", NA),
      SA      = c("SA", "nonSA", NA)
    ),

    # time × SA and time × type interactions
    expand_grid(diff = "85 - 30", diff2 = "SA - nonSA",         timebin = NA, actype = c("mental", "physical", NA), SA = NA),
    expand_grid(diff = "85 - 30", diff2 = "physical - mental",  timebin = NA, actype = NA,                         SA = c("SA", "nonSA", NA))

  ) %>%
    do.call(rbind.data.frame, .) %>%
    mutate(
      SA   = case_when(
        diff == "physical_SA - mental_SA"       ~ "SA",
        diff == "physical_nonSA - mental_nonSA" ~ "nonSA",
        .default = SA
      ),
      diff = gsub("_SA|_nonSA", "", diff)
    )

  ## ---- Compute contrasts ----
  contrasts <- lapply(
    set_names(x = c("log", "logexp")),
    function(s) with(
      post_comp,
      sapply(
        X   = seq_len(nrow(post_comp)),
        FUN = function(i) {

          term1 <- strsplit(diff[i], " - ")[[1]][1]
          term2 <- strsplit(diff[i], " - ")[[1]][2]

          # identify which dimension is being compared
          if (term1 %in% na.omit(unique(actype))) {
            col1 <- paste(SA[i], term1, timebin[i], sep = "_")
            col2 <- paste(SA[i], term2, timebin[i], sep = "_")
          } else if (term1 %in% na.omit(unique(SA))) {
            col1 <- paste(term1, actype[i], timebin[i], sep = "_")
            col2 <- paste(term2, actype[i], timebin[i], sep = "_")
          } else if (term1 %in% na.omit(unique(timebin))) {
            col1 <- paste(SA[i], actype[i], term1, sep = "_")
            col2 <- paste(SA[i], actype[i], term2, sep = "_")
          }

          dif <- if (s == "logexp") exp(ppred[, col1]) - exp(ppred[, col2])
          else                   ppred[, col1]  -     ppred[, col2]

          # two-way interaction of contrasts
          if (!is.na(diff2[i])) {
            term21 <- strsplit(diff2[i], " - ")[[1]][1]
            term22 <- strsplit(diff2[i], " - ")[[1]][2]

            if (term1 %in% na.omit(unique(actype))) {
              col11 <- paste(term21, term1, timebin[i], sep = "_")
              col12 <- paste(term21, term2, timebin[i], sep = "_")
              col21 <- paste(term22, term1, timebin[i], sep = "_")
              col22 <- paste(term22, term2, timebin[i], sep = "_")
            } else if (term1 %in% na.omit(unique(timebin))) {
              if (term21 %in% na.omit(unique(actype))) {
                col11 <- paste(SA[i], term21, term1, sep = "_")
                col12 <- paste(SA[i], term21, term2, sep = "_")
                col21 <- paste(SA[i], term22, term1, sep = "_")
                col22 <- paste(SA[i], term22, term2, sep = "_")
              } else if (term21 %in% na.omit(unique(SA))) {
                col11 <- paste(term21, actype[i], term1, sep = "_")
                col12 <- paste(term21, actype[i], term2, sep = "_")
                col21 <- paste(term22, actype[i], term1, sep = "_")
                col22 <- paste(term22, actype[i], term2, sep = "_")
              }
            }

            dif <- if (s == "logexp") {
              (exp(ppred[, col11]) - exp(ppred[, col12])) - (exp(ppred[, col21]) - exp(ppred[, col22]))
            } else {
              (ppred[, col11] - ppred[, col12]) - (ppred[, col21] - ppred[, col22])
            }
          }

          est <- paste0(rprint(median(dif), 2), " ", ciprint(unlist(eti(dif))[-1]))
          pd  <- paste0(rprint(100 * c(p_direction(dif))$pd, 2), "%")
          return(c(Estimate = est, pd = pd))
        }
      )
    ) %>%
      t() %>%
      cbind.data.frame(post_comp, .) %>%
      mutate(scale = s, .before = Estimate)
  ) %>%
    do.call(rbind.data.frame, .) %>%
    pivot_wider(names_from = scale, values_from = c(Estimate, pd)) %>%
    relocate(pd_log, .after = Estimate_log)

  ## ---- Return requested object ----
  return(get(output))
}


#' Plot posterior conditional means as interaction plots
#'
#' @description
#' Produces a set of \pkg{ggplot2} interaction plots from the posterior
#' expectations data frame returned by [compute_posterior_expectations()].
#' For each scale (`"log"` and back-transformed `"explog"`), two complementary
#' panels are drawn:
#' \enumerate{
#'   \item Activity type on colour, SuperAging status on facet rows.
#'   \item SuperAging status on colour, activity type on facet rows.
#' }
#' The two panels are combined with \pkg{patchwork}. Optionally, an
#' informative title and subtitle are added, and the plot is saved as a JPEG.
#'
#' @param posterior_expect A data frame as returned by
#'   [compute_posterior_expectations()] with `output = "expectations"`.
#' @param text Logical. If `TRUE`, adds a title and subtitle to the combined
#'   plot via [patchwork::plot_annotation()]. Defaults to `FALSE`.
#' @param save Logical. If `TRUE` (default), saves the combined plot to
#'   `"figures/conditional_means_<scale>_scale.jpg"` at 300 dpi.
#'
#' @return A named list with elements `"log"` and `"back-transformed raw"`,
#'   each a `patchwork` object (or `NULL` if plotting failed). Called
#'   primarily for its side effect of saving figures when `save = TRUE`.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_linerange geom_line
#'   scale_colour_manual labs facet_wrap theme_bw theme element_text
#'   ggsave last_plot
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom here here
#'
#' @seealso [compute_posterior_expectations()]
#'
#' @export
draw_interaction_plots <- function(posterior_expect, text = FALSE, save = TRUE) {

  # build individual panels for each scale
  plts <- lapply(
    set_names(x = c("log", "explog"), nm = c("log", "back-transformed raw")),
    function(i) {

      data <- posterior_expect %>%
        filter(!is.na(Activity_type), !is.na(SA), scale == i) %>%
        mutate(
          SuperAging = factor(
            x      = if_else(SA == "SA", 1, 0),
            levels = 0:1,
            labels = c("Yes", "No")
          ),
          `Activity type` = factor(
            x      = case_when(
              Activity_type == "physical" ~ 1,
              Activity_type == "mental"   ~ 2
            ),
            levels = 1:2,
            labels = c("Physical", "Mental")
          )
        )

      list(
        # panel A: colour = activity type, facet = SA
        data %>%
          ggplot() +
          aes(x = Time_bin, y = Estimate, ymin = ETI_low, ymax = ETI_high,
              colour = `Activity type`, group = `Activity type`) +
          geom_point(position = position_dodge(width = .5), size = 4) +
          geom_linerange(position = position_dodge(width = .5), linewidth = 1.5) +
          geom_line(linetype = "dashed", linewidth = .8) +
          scale_colour_manual(values = c("#56B4E9", "#E69F00")) +
          labs(
            y = ifelse(i == "log", "log(Intensity)", "Intensity (1\u20135 points)"),
            x = "Time bin (years)"
          ) +
          facet_wrap(~ SA, nrow = 2,
                     labeller = as_labeller(c(SA = "SuperAging", nonSA = "non-SuperAging"))) +
          theme_bw(base_size = 12) +
          theme(legend.position = "bottom"),

        # panel B: colour = SA, facet = activity type
        data %>%
          ggplot() +
          aes(x = Time_bin, y = Estimate, ymin = ETI_low, ymax = ETI_high,
              colour = SuperAging, group = SuperAging) +
          geom_point(position = position_dodge(width = .5), size = 4) +
          geom_linerange(position = position_dodge(width = .5), linewidth = 1.5) +
          geom_line(linetype = "dashed", linewidth = .8) +
          scale_colour_manual(values = c("#CC79A7", "#999999")) +
          labs(
            y = ifelse(i == "log", "log(Intensity)", "Intensity (1\u20135 points)"),
            x = "Time bin (years)"
          ) +
          facet_wrap(~ `Activity type`, nrow = 2) +
          theme_bw(base_size = 12) +
          theme(legend.position = "bottom")
      )
    }
  )

  # combine panels and optionally annotate / save
  figs <- lapply(
    set_names(x = names(plts)),
    function(i) {

      output <- with(plts, get(i)[[1]] | get(i)[2]) +
        plot_layout(axis_titles = "collect")

      if (text) {
        output <- output +
          plot_annotation(
            title    = "Three-way interaction between time bin, SA, and activity type",
            subtitle = paste0("Left and right sets of panels are representations of the same interaction on a ", i, " scale"),
            theme    = theme(
              plot.title    = element_text(hjust = .5, face = "bold"),
              plot.subtitle = element_text(hjust = .5)
            )
          )
      }

      new_folder("figures")
      fn <- paste0("conditional_means_", gsub("-| ", "_", i), "_scale.jpg")

      if (save) {
        ggsave(
          plot     = last_plot(),
          filename = here("figures", fn),
          dpi      = 300,
          width    = 12.6,
          height   = 8.31
        )
      }
    }
  )

  return(figs)
}
