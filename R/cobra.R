# ── COBRA-A mental activity ───────────────────────────────────────────────────

#' Read and pre-process COBRA-A mental activity data
#'
#' @description
#' Reads the full `info.rds` participant metadata file and the
#' `cobra_a.rds` activity-level records, derives a SuperAger factor
#' (`saf`), restricts the activity set to activities with complete data for
#' every participant, and builds a per-participant summary data frame (`data`)
#' together with the filtered long-format activity records (`cad`).
#'
#' Three per-participant summary variables are computed:
#' \describe{
#'   \item{`n_activity`}{Number of activities reported at least once per
#'     week (`times_per_week > 0`).}
#'   \item{`total_hours_per_week`}{Total hours across all activities.}
#'   \item{`total_times_per_week`}{Total weekly occasions across all
#'     activities.}
#' }
#'
#' @param cobra_a A data frame of COBRA-A long-format activity records with at
#'   minimum the columns `id` (character), `activity` (character),
#'   `times_per_week` (numeric), and `hours_per_week` (numeric). Typically
#'   produced by `readRDS("data-raw/cobra_a.rds")`.
#' @param info A data frame of participant metadata as returned by
#'   `readRDS("data-raw/info.rds")`, containing at minimum the columns `id`
#'   (character), `sa` (logical SuperAging flag), `age` (numeric),
#'   `mmse` (numeric), and `edu2` (integer education level, 1 or 2).
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{`data`}{A per-participant data frame (`dd`) with columns `id`,
#'       `n_activity`, `total_hours_per_week`, `total_times_per_week`, `sa`,
#'       `saf`, `age`, `mmse`, and `edu2`. Participants with missing `sa`
#'       are excluded.}
#'     \item{`cad`}{The long-format activity records restricted to the
#'       complete-data activities in `activities`.}
#'     \item{`activities`}{A character vector of activity names present for
#'       every participant (used to define the analysis set).}
#'   }
#'
#' @importFrom dplyr count filter pull mutate group_by summarise left_join
#'   select join_by
#'
#' @seealso [fit_cobra_regressions()], [run_cobra_activity_tests()]
#'
#' @export
build_cobra_analysis_data <- function(cobra_a, info) {

  nid <- length(unique(cobra_a$id))

  # restrict to activities with data for every participant
  cobra_activity <- cobra_a |>
    count(activity) |>
    filter(n == nid) |>
    pull(activity)

  cad <- cobra_a |> filter(activity %in% cobra_activity)

  info_aug <- info |>
    mutate(
      saf = factor(
        sa,
        levels = c(TRUE, FALSE),
        labels = c("Superager", "non-Superager")
      )
    )

  dd <- cad |>
    group_by(id) |>
    summarise(
      n_activity           = sum(times_per_week > 0),
      total_hours_per_week = sum(hours_per_week),
      total_times_per_week = sum(times_per_week),
      .groups = "drop"
    ) |>
    left_join(
      info_aug |> select(id, sa, saf, age, mmse, edu2),
      by = join_by(id)
    ) |>
    filter(!is.na(sa))

  list(data = dd, cad = cad, activities = cobra_activity)
}


#' Fit regression models for COBRA-A mental activity outcomes
#'
#' @description
#' Fits three OLS regression models predicting, respectively, the number of
#' activities (`n_activity`), total hours per week (`total_hours_per_week`),
#' and total weekly occasions (`total_times_per_week`) from SuperAger status
#' (`saf`), age, and education level (`edu2`).
#'
#' The coefficient matrix `cbind(coef(m), confint(m))` is precomputed for
#' each model and returned alongside the fitted objects so that inline
#' reporting (e.g. via [rndci()]) does not need to re-call `confint()`.
#'
#' @param data The `$data` element of the list returned by
#'   [build_cobra_analysis_data()]: a per-participant data frame with columns
#'   `n_activity`, `total_hours_per_week`, `total_times_per_week`, `saf`,
#'   `age`, and `edu2`.
#'
#' @return A named list with six elements, two per outcome:
#'   \describe{
#'     \item{`n_activity`, `hours`, `times`}{Fitted `lm` objects.}
#'     \item{`coef_n_activity`, `coef_hours`, `coef_times`}{Numeric matrices
#'       with one row per predictor and three columns: `Estimate`, `2.5 %`,
#'       `97.5 %`.}
#'   }
#'
#' @seealso [build_cobra_analysis_data()], [rndci()]
#'
#' @export
fit_cobra_regressions <- function(data) {

  m1 <- lm(n_activity           ~ saf + age + edu2, data = data)
  m2 <- lm(total_hours_per_week ~ saf + age + edu2, data = data)
  m3 <- lm(total_times_per_week ~ saf + age + edu2, data = data)

  list(
    n_activity       = m1,
    hours            = m2,
    times            = m3,
    coef_n_activity  = cbind(Estimate = coef(m1), confint(m1)),
    coef_hours       = cbind(Estimate = coef(m2), confint(m2)),
    coef_times       = cbind(Estimate = coef(m3), confint(m3))
  )
}


#' Run an independent-samples t-test for a single COBRA-A activity
#'
#' @description
#' Filters the long-format COBRA-A activity records to a single activity,
#' aggregates total hours per week per participant, joins SuperAger status,
#' and passes the two group vectors to [apa::t_test()].
#'
#' @param activity_name A length-one character string naming the activity to
#'   test (must be present in `cobra_analysis$cad$activity`).
#' @param cobra_analysis The named list returned by
#'   [build_cobra_analysis_data()], containing elements `cad` and `data`.
#'
#' @return An object of class `htest` as returned by [apa::t_test()], with an
#'   additional `cohens_d` attribute produced by [apa::cohens_d()]. The
#'   object can be formatted with [apa::t_apa()].
#'
#' @importFrom dplyr filter group_by summarise left_join select join_by
#'
#' @seealso [run_cobra_activity_tests()], [apa::t_test()], [apa::t_apa()]
#'
#' @export
test_cobra_activity <- function(activity_name, cobra_analysis) {

  d1 <- cobra_analysis$cad |>
    filter(activity == activity_name) |>
    group_by(id) |>
    summarise(total_hours_per_week = sum(hours_per_week), .groups = "drop") |>
    left_join(
      cobra_analysis$data |> select(id, sa),
      by = join_by(id)
    ) |>
    filter(!is.na(sa))

  apa::t_test(
    d1$total_hours_per_week[d1$sa],
    d1$total_hours_per_week[!d1$sa]
  )
}


#' Run t-tests for all COBRA-A activities
#'
#' @description
#' Applies [test_cobra_activity()] to every activity in
#' `cobra_analysis$activities` and returns the results as a named list.
#' The list is suitable for passing to [make_cobra_activity_table()] and for
#' direct inline reporting in a Quarto / R Markdown document.
#'
#' @param cobra_analysis The named list returned by
#'   [build_cobra_analysis_data()].
#'
#' @return A named list of `htest` objects (one per activity), where names
#'   correspond to `cobra_analysis$activities`.
#'
#' @seealso [test_cobra_activity()], [make_cobra_activity_table()]
#'
#' @export
run_cobra_activity_tests <- function(cobra_analysis) {
  setNames(
    lapply(cobra_analysis$activities, test_cobra_activity,
           cobra_analysis = cobra_analysis),
    cobra_analysis$activities
  )
}


#' Build a per-activity summary table for COBRA-A mental activities
#'
#' @description
#' For each activity, computes mean and SD of `hours_per_week` separately for
#' SuperAgers and non-SuperAgers, merges in the pre-computed t-test results
#' (from [run_cobra_activity_tests()]), and formats the result as a \pkg{gt}
#' table sorted by decreasing absolute Cohen's *d*.
#'
#' Columns are arranged as: Activity | SA Mean (SD) | non-SA Mean (SD) |
#' *d* | *p*.
#'
#' @param cobra_analysis The named list returned by
#'   [build_cobra_analysis_data()].
#' @param activity_tests A named list of `htest` objects as returned by
#'   [run_cobra_activity_tests()].
#'
#' @return A `gt_tbl` object.
#'
#' @importFrom dplyr left_join filter group_by summarise mutate arrange
#'   select join_by
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map map_dbl
#' @importFrom stringr str_replace_all
#' @importFrom gt gt fmt_number cols_label tab_spanner tab_footnote md
#' @importFrom gtExtras cols_merge_n_pct
#'
#' @seealso [run_cobra_activity_tests()], [build_cobra_analysis_data()]
#'
#' @export
make_cobra_activity_table <- function(cobra_analysis, activity_tests) {

  tab <- cobra_analysis$cad |>
    left_join(
      cobra_analysis$data |> select(id, sa, saf, age, mmse, edu2),
      by = join_by(id)
    ) |>
    filter(!is.na(saf)) |>
    group_by(activity, saf) |>
    summarise(
      m  = mean(hours_per_week),
      sd = sd(hours_per_week),
      n  = n(),
      .groups = "drop"
    ) |>
    pivot_wider(
      id_cols     = "activity",
      names_from  = "saf",
      values_from = c("m", "sd", "n")
    ) |>
    mutate(
      tt = map(activity, ~ activity_tests[[.x]]),
      d  = map_dbl(tt, apa::cohens_d),
      p  = map_dbl(tt, \(x) x$p.value)
    ) |>
    arrange(-abs(d))

  tab |>
    select(activity, `m_Superager`, `n_Superager`,
           `m_non-Superager`, `n_non-Superager`, d, p) |>
    mutate(activity = str_replace_all(activity, "_", " ")) |>
    #gt() |>
    gt_apa(title = "Average time spent in each activity per week") |>
    fmt_number(columns = 2:5, decimals = 1) |>
    fmt_number(columns = 6,   decimals = 3) |>
    fmt_number(columns = 7,   decimals = 3) |>
    cols_merge_n_pct(col_n = 2, col_pct = 3) |>
    cols_merge_n_pct(col_n = 4, col_pct = 5) |>
    tab_spanner(label = "Superagers",     columns = 2) |>
    tab_spanner(label = "non-Superagers", columns = 4) |>
    cols_label(
      .list = list(
        "m_Superager"     = "Mean (SD)",
        "m_non-Superager" = "Mean (SD)",
        "d"               = md("*d*"),
        "p"               = md("*p*"),
        "activity"        = "Activity"
      )
    ) |>
    tab_footnote(
      "Activities are sorted by decreasing absolute Cohen's d."
    )
}


# ── Physical activity (WHO criteria) ──────────────────────────────────────────

#' Analyse WHO physical activity criteria in relation to SuperAging
#'
#' @description
#' Joins the COBRA-B / WHO physical activity data with the participant
#' metadata, and runs a chi-square test of independence between SuperAger
#' status and meeting WHO physical activity guidelines.
#'
#' @param info A data frame as returned by `readRDS("data-raw/info.rds")`,
#'   containing at minimum the columns `id` (character), `sa` (logical), and
#'   `saf` (factor, derived by [build_cobra_analysis_data()] or elsewhere).
#'   If `saf` is absent it is computed from `sa`.
#' @param cobra_b A data frame as returned by
#'   `readRDS("data-raw/cobra_b_who.rds")`, containing at minimum the columns
#'   `id` (character), `who_pa` (integer 0/1, whether the participant met WHO
#'   criteria), and `who_paf` (factor version of `who_pa`).
#'
#' @return A named list with four elements:
#'   \describe{
#'     \item{`data`}{The joined data frame (`db`) used for all tests.}
#'     \item{`n_who_pa`}{Integer. Number of participants meeting WHO criteria.}
#'     \item{`pct_who_pa`}{Numeric. Percentage meeting WHO criteria (0–100).}
#'     \item{`chisq`}{An `htest` object from [stats::chisq.test()],
#'       comparing `saf` × `who_paf`.}
#'   }
#'
#' @importFrom dplyr mutate left_join join_by
#'
#' @seealso [build_cobra_analysis_data()]
#'
#' @export
compare_physical_activity <- function(info, cobra_b) {

  db <- info |>
    mutate(
      saf = if ("saf" %in% names(info)) saf else
        factor(sa, levels = c(TRUE, FALSE), labels = c("Superager", "non-Superager"))
    ) |>
    left_join(cobra_b, by = join_by(id))

  list(
    data       = db,
    n_who_pa   = sum(db$who_pa,  na.rm = TRUE),
    pct_who_pa = mean(db$who_pa, na.rm = TRUE) * 100,
    chisq      = chisq.test(table(db$saf, db$who_paf))
  )
}


# ── Social activity (VLS-AQL) ─────────────────────────────────────────────────

#' Compare social activity (VLS-AQL) between SuperAging groups
#'
#' @description
#' Joins VLS Activity of Later Life Questionnaire (VLS-AQL) scores with
#' participant metadata and tests for differences between SuperAgers and
#' non-SuperAgers in both the private and public social activity subscales
#' using independent-samples *t*-tests (via [apa::t_test()]).
#'
#' @param info A data frame as returned by `readRDS("data-raw/info.rds")`,
#'   containing at minimum the columns `id` (character) and `sa` (logical).
#' @param vls A data frame as returned by `readRDS("data-raw/vls.rds")`,
#'   where the first three columns are `id` (character),
#'   `vls_alq_private` (numeric), and `vls_alq_public` (numeric).
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{`data`}{The joined data frame (`dv`) used for all tests.}
#'     \item{`means`}{A named numeric vector with four elements:
#'       `sa_private`, `nonsa_private`, `sa_public`, `nonsa_public`.}
#'     \item{`private_ttest`}{An `apa::t_test` `htest` object for the
#'       private subscale.}
#'     \item{`public_ttest`}{An `apa::t_test` `htest` object for the
#'       public subscale.}
#'   }
#'
#' @importFrom dplyr left_join select join_by
#'
#' @seealso [apa::t_test()], [apa::t_apa()]
#'
#' @export
compare_social_activity <- function(info, vls) {

  dv <- info |>
    left_join(vls |> select(1:3), by = join_by(id))

  priv <- dv$vls_alq_private
  pub  <- dv$vls_alq_public
  sa   <- dv$sa

  list(
    data = dv,
    means = c(
      sa_private    = mean(priv[sa],   na.rm = TRUE),
      nonsa_private = mean(priv[!sa],  na.rm = TRUE),
      sa_public     = mean(pub[sa],    na.rm = TRUE),
      nonsa_public  = mean(pub[!sa],   na.rm = TRUE)
    ),
    private_ttest = apa::t_test(priv[sa],  priv[!sa]),
    public_ttest  = apa::t_test(pub[sa],   pub[!sa])
  )
}
