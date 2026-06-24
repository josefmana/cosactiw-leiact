#' Count participants by SuperAging status
#'
#' @description
#' Returns the number of participants in each SuperAging (`"SA"`) and
#' non-SuperAging (`"nonSA"`) group.
#'
#' @param data A data frame containing at minimum a column `SA` with values
#'   `"SA"` and `"nonSA"`.
#'
#' @return A named integer vector of length two with names `"SA"` and
#'   `"nonSA"`, each element giving the row count for that group.
#'
#' @examples
#' \dontrun{
#' count_subjects(wide_data)
#' }
#'
#' @export
count_subjects <- function(data) {
  sapply(
    c("SA", "nonSA"),
    \(i) nrow(subset(data, SA == i))
  )
}


#' Compare demographic variables between SuperAging groups
#'
#' @description
#' Computes means and standard deviations (via [cenvar()]) of `Age_years` and
#' `Education_years` separately for SuperAgers and non-SuperAgers, then
#' appends the results of independent-samples *t*-tests (statistic, degrees of
#' freedom, and *p*-value) for each variable.
#'
#' @param data A data frame containing at minimum the columns `SA`
#'   (`"SA"` / `"nonSA"`), `Age_years` (numeric), and `Education_years`
#'   (numeric).
#'
#' @return A data frame with one row per demographic variable and columns:
#'   \describe{
#'     \item{`SA`, `nonSA`}{Formatted mean ± SD strings for each group.}
#'     \item{`t`}{Formatted *t*-statistic (3 d.p.).}
#'     \item{`df`}{Formatted degrees of freedom (2 d.p.).}
#'     \item{`p`}{Formatted *p*-value string.}
#'   }
#'
#' @importFrom dplyr group_by summarise across mutate
#' @importFrom tibble column_to_rownames
#'
#' @seealso [cenvar()], [rprint()], [pprint()]
#'
#' @export
compare_demography <- function(data) {
  data %>%
    group_by(SA) %>%
    summarise(across(all_of(c("Age_years", "Education_years")), cenvar)) %>%
    column_to_rownames("SA") %>%
    t() %>%
    as.data.frame() %>%
    mutate(
      t  = unlist(
        sapply(rownames(.), function(y)
          rprint(t.test(as.formula(paste0(y, " ~ SA")), data = data)$statistic, d = 3)),
        use.names = FALSE
      ),
      df = unlist(
        sapply(rownames(.), function(y)
          rprint(t.test(as.formula(paste0(y, " ~ SA")), data = data)$parameter, d = 2)),
        use.names = FALSE
      ),
      p  = unlist(
        sapply(rownames(.), function(y)
          pprint(t.test(as.formula(paste0(y, " ~ SA")), data = data)$p.value, dec = 3, text = TRUE)),
        use.names = FALSE
      )
    )
}


#' Build a regression-type specification table
#'
#' @description
#' Constructs a data frame enumerating all outcome variables to be modelled
#' (the total activity count, each activity type, and each activity category)
#' together with the regression family to use for each (`"lm"` or `"glm"`).
#' This specification table is consumed by [count_regressions()].
#'
#' @param input A named list as returned by [get_data()], containing at
#'   minimum the element `map` (activity mapping table with columns `type` and
#'   `category`).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{`y`}{Character. Name of the outcome variable.}
#'     \item{`type`}{Character. Regression type, currently always `"lm"`.}
#'   }
#'
#' @seealso [count_regressions()]
#'
#' @export
specify_regression <- function(input) {
  data.frame(
    y = with(input, c("activities", unique(map$type), unique(map$category)))
  ) %>%
    mutate(
      type = case_when(
        y %in% c("activities") ~ "lm",
        .default = "lm"
      )
    )
}


#' Fit linear (or Poisson GLM) regressions for activity counts
#'
#' @description
#' Fits a separate regression model for each outcome variable listed in
#' `specs`, predicting activity counts from SuperAging status, age, and a
#' binary education factor. Predictors are recoded internally to match the
#' conventions used in the VLS and COBRA analyses:
#' \itemize{
#'   \item `saf`: `factor` — `"Superager"` vs `"non-Superager"`.
#'   \item `edu2`: `factor` — education level dichotomised at level 2.
#'   \item `age`: numeric age in years.
#' }
#' The regression family is read from `specs$type`; `"lm"` calls [stats::lm()]
#' and `"glm"` fits a Poisson log-linear model via [stats::glm()].
#'
#' @param data A wide-format data frame as returned by [pivot_data()],
#'   containing the outcome columns named in `specs$y` as well as `SA`,
#'   `Age_years`, and `Education_level`.
#' @param specs A data frame as returned by [specify_regression()] with
#'   columns `y` (outcome name) and `type` (`"lm"` or `"glm"`).
#'
#' @return A named list of fitted model objects (one per row of `.specs`),
#'   where names correspond to the outcome variable names in `.specs$y`.
#'
#' @importFrom dplyr mutate if_else
#'
#' @seealso [specify_regression()], [diagnose_regressions()],
#'   [extract_parameters()]
#'
#' @export
count_regressions <- function(data, specs) {

  # recode variables to match VLS / COBRA conventions
  d0 <- data %>%
    mutate(
      saf = factor(
        if_else(SA == "SA", TRUE, FALSE),
        levels = c(TRUE, FALSE),
        labels = c("Superager", "non-Superager")
      ),
      edu2 = factor(
        if_else(as.numeric(Education_level) > 2, 2, 1),
        levels = 1:2
      ),
      age = Age_years
    )

  # fit models
  fits <- with(
    specs,
    lapply(
      set_names(y),
      function(i)
        if (type[y == i] == "glm") {
          glm(
            formula = as.formula(paste0(i, " ~ saf + age + edu2")),
            data    = d0,
            family  = poisson(link = "log")
          )
        } else {
          do.call(
            what = type[y == i],
            args = list(
              formula = as.formula(paste0(i, " ~ saf + age + edu2")),
              data    = d0
            )
          )
        }
    )
  )

  return(fits)
}


#' Run model-diagnostic checks on a list of regression fits
#'
#' @description
#' Applies [performance::check_model()] to each element of a named list of
#' fitted regression objects, producing a set of diagnostic plots for visual
#' inspection of model assumptions (linearity, homoscedasticity, normality of
#' residuals, influential observations, etc.).
#'
#' @param fits A named list of fitted model objects as returned by
#'   [count_regressions()].
#'
#' @return A named list of `check_model` plot objects (one per model), with
#'   the same names as `fits`.
#'
#' @importFrom performance check_model
#'
#' @seealso [count_regressions()], [performance::check_model()]
#'
#' @export
diagnose_regressions <- function(fits) {
  lapply(
    set_names(names(fits)),
    function(y) check_model(x = fits[[y]])
  )
}


#' Summarise activity counts by SuperAging status
#'
#' @description
#' Computes formatted mean ± SD strings (via [cenvar()]) for each activity
#' outcome variable (total count, type totals, category totals) separately for
#' SuperAgers and non-SuperAgers. Outcome variables with a zero sum in a
#' subgroup are represented as `"-"`.
#'
#' @param input A named list as returned by [get_data()], used to determine
#'   the outcome variable names from `map$type` and `map$category`.
#' @param data A wide-format data frame as returned by [pivot_data()],
#'   containing the outcome columns and the `SA` column.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{`type`}{Character. Broad activity domain: `"Both"`, `"Mental"`,
#'       or `"Physical"`.}
#'     \item{`category`}{Character. Activity category name, or `"all"` for
#'       domain totals.}
#'     \item{`nonSA`}{Formatted mean ± SD for non-SuperAgers (or `"-"`).}
#'     \item{`SA`}{Formatted mean ± SD for SuperAgers (or `"-"`).}
#'   }
#'   Rows are arranged by `type`.
#'
#' @seealso [cenvar()], [add_comparisons()]
#'
#' @export
count_activities <- function(input, data) {

  with(input, lapply(
    c("activities", unique(map$type), unique(map$category)),
    function(y) sapply(
      c("nonSA", "SA"),
      function(i) if_else(
        condition = sum(subset(data, SA == i)[, y]) == 0,
        true  = "-",
        false = paste0(cenvar(y = subset(data, SA == i)[, y], dec = 1, cen = "mean", var = "sd", sep = " ("), ")")
      )
    ) %>%
      t() %>%
      as.data.frame() %>%
      mutate(category = y, .before = 1)
  ) %>%
    do.call(rbind.data.frame, .) %>%
    mutate(
      type = sapply(
        seq_len(nrow(.)),
        function(i) case_when(
          category[i] == "activities"                                                                         ~ "all",
          category[i] == "mental"   | category[i] %in% unique(subset(input$map, type == "mental")$category)   ~ "mental",
          category[i] == "physical" | category[i] %in% unique(subset(input$map, type == "physical")$category) ~ "physical"
        )
      ),
      category = case_when(
        category %in% c("activities", "mental", "physical") ~ "all",
        .default = category
      ),
      .before = 1
    ) %>%
    arrange(type))
}


#' Append group-comparison statistics to an activity-count table (main text)
#'
#' @description
#' For each row of a descriptive table produced by [count_activities()],
#' determines the corresponding outcome variable, runs an independent-samples
#' *t*-test between SuperAging groups, and appends Cohen's *d*, the raw
#' mean difference with its 95% CI, and the *p*-value as new columns.
#'
#' The difference is in the SA − non-SA direction: positive values indicate
#' SuperAgers report higher counts.  The formula interface is used throughout,
#' so group ordering follows R's alphabetical convention (\"SA\" precedes
#' \"nonSA\", giving estimate[[1]] = mean(SA)).
#'
#' @param table A data frame as returned by [count_activities()], with
#'   columns `type` and `category`.
#' @param data A wide-format data frame as returned by [pivot_data()],
#'   containing the outcome columns and the `SA` column.
#'
#' @return A data frame combining `table` with the following additional
#'   columns:
#'   \describe{
#'     \item{`cohens_d`}{Formatted Cohen's *d*.}
#'     \item{`diff_ci`}{Raw SA − non-SA difference with 95% CI, formatted
#'       as `\"estimate [lower, upper]\"`.)
#'     \item{`p_ttest`}{Formatted *t*-test *p*-value.}
#'   }
#'   `NaN` and `NA` values are replaced with `"-"`.
#'
#' @importFrom effsize cohen.d
#'
#' @seealso [count_activities()], [add_comparisons_supp()],
#'   [add_regression_parameters()]
#'
#' @export
add_comparisons <- function(table, data) {

  sapply(
    seq_len(nrow(table)),
    function(i) {

      # determine the outcome variable for this row
      v <- with(table, case_when(
        category[i] == "all" & type[i] == "all" ~ "activities",
        category[i] == "all"                    ~ type[i],
        .default                                = category[i]
      ))

      # parametric test
      t <- t.test(as.formula(paste0(v, " ~ SA")), data = data)
      d <- cohen.d(as.formula(paste0(v, " ~ SA")), data = data)

      # SA - nonSA difference: "SA" < "nonSA" alphabetically so
      # estimate[[1]] = mean(SA), estimate[[2]] = mean(nonSA)
      diff_ci <- sprintf(
        "%s [%s, %s]",
        rprint(t$estimate[[1]] - t$estimate[[2]], 1),
        rprint(t$conf.int[[1]], 1),
        rprint(t$conf.int[[2]], 1)
      )

      return(c(
        diff_ci  = diff_ci,
        cohens_d = rprint(d$estimate, 3),
        t        = rprint(t$statistic, 2), # keep in for in-text shenaningans
        df       = rprint(t$parameter, 2), # keep in for in-text shenaningans
        p        = zerolead(t$p.value, 3)
      ))
    }
  ) %>%
    t() %>%
    as.data.frame() %>%
    cbind.data.frame(table, .) %>%
    mutate(
      across(everything(), ~ case_when(
        grepl("NaN", .x) ~ "-",
        is.na(.x)        ~ "-",
        .default = .x
      )),
      category = if_else(category == "all", "-", category),
      type = case_when(
        type == "all"      ~ "Both",
        type == "mental"   ~ "Mental",
        type == "physical" ~ "Physical"
      )
    )
}


#' Append non-parametric comparison statistics to an activity-count table
#' (appendix)
#'
#' @description
#' Non-parametric counterpart of [add_comparisons()] intended for the
#' supplementary appendix.  For each row of a descriptive table produced by
#' [count_activities()], computes the median and range (min–max) for each
#' SuperAging group, runs a Wilcoxon rank-sum test, and computes the Vargha–
#' Delaney *A* effect size.
#'
#' VD.A is computed in the SA > non-SA direction (values > 0.5 indicate
#' SuperAgers tend to report higher counts).
#'
#' @param table A data frame as returned by [count_activities()], with
#'   columns `type` and `category`.  Only the structural columns are carried
#'   forward; descriptive statistics are recomputed from `data`.
#' @param data A wide-format data frame as returned by [pivot_data()],
#'   containing the outcome columns and the `SA` column.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{`type`}{Activity domain, recoded to `"Both"`, `"Mental"`, or
#'       `"Physical"`.}
#'     \item{`category`}{Activity category, or `"-"` for domain totals.}
#'     \item{`SA_med_range`}{Formatted median (min–max) for SuperAgers.}
#'     \item{`nonSA_med_range`}{Formatted median (min–max) for non-SuperAgers.}
#'     \item{`VDA`}{Formatted Vargha–Delaney *A*.}
#'     \item{`W`}{Formatted Wilcoxon *W* statistic.}
#'     \item{`p_mannwhitney`}{Formatted *p*-value.}
#'   }
#'   `NaN` and `NA` values are replaced with `"-"`.
#'
#' @importFrom effsize VD.A
#'
#' @seealso [count_activities()], [add_comparisons()]
#'
#' @export
add_comparisons_supp <- function(table, data) {

  fmt_med_range <- function(x) {
    sprintf(
      "%s (%s\u2013%s)",
      rprint(median(x, na.rm = TRUE), 1),
      rprint(min(x,    na.rm = TRUE), 0),
      rprint(max(x,    na.rm = TRUE), 0)
    )
  }

  sapply(
    seq_len(nrow(table)),
    function(i) {

      # determine the outcome variable for this row
      v <- with(table, case_when(
        category[i] == "all" & type[i] == "all" ~ "activities",
        category[i] == "all"                    ~ type[i],
        .default                                = category[i]
      ))

      sa_vals    <- data[[v]][data$SA == "SA"]
      nonsa_vals <- data[[v]][data$SA == "nonSA"]

      # guard against degenerate (all-zero) rows present in some categories
      if (sum(sa_vals, na.rm = TRUE) == 0 && sum(nonsa_vals, na.rm = TRUE) == 0) {
        return(c(
          SA_med_range    = "-",
          nonSA_med_range = "-",
          VDA             = "-",
          W               = "-",
          p               = "-"
        ))
      }

      wt  <- wilcox.test(sa_vals, nonsa_vals)
      vda <- VD.A(sa_vals, nonsa_vals)

      return(c(
        SA_med_range    = fmt_med_range(sa_vals),
        nonSA_med_range = fmt_med_range(nonsa_vals),
        VDA             = rprint(vda$estimate, 3),
        W               = rprint(wt$statistic, 1),
        p               = pprint(wt$p.value, 3, text = FALSE)
      ))
    }
  ) %>%
    t() %>%
    as.data.frame() %>%
    cbind.data.frame(
      table %>% select(type, category),  # structural columns only
      .
    ) %>%
    mutate(
      across(everything(), ~ case_when(
        grepl("NaN", .x) ~ "-",
        is.na(.x)        ~ "-",
        .default = .x
      )),
      category = if_else(category == "all", "-", category),
      type = case_when(
        type == "all"      ~ "Both",
        type == "mental"   ~ "Mental",
        type == "physical" ~ "Physical"
      )
    )
}



#'
#' @description
#' For each model in `fits`, extracts the coefficient estimates, 95%
#' confidence intervals, *t*-values, and *p*-values, then pivots the result
#' to a wide format suitable for joining with a descriptive table. Each
#' predictor's statistics are stored in a group of four columns named
#' `b_<predictor>`, `CI_<predictor>`, `t_<predictor>`, and `p_<predictor>`.
#'
#' @param fits A named list of `lm` fit objects as returned by
#'   [count_regressions()].
#' @param input A named list as returned by [get_data()], used to classify
#'   outcome variables by `type` and `category` via `input$map`.
#'
#' @return A data frame with columns `type` and `category` (for joining) plus
#'   four formatted columns per predictor in each model. The predictor
#'   `(Intercept)` is included.
#'
#' @importFrom stats coefficients confint.lm summary
#' @importFrom tidyr pivot_wider
#' @importFrom tibble rownames_to_column
#'
#' @seealso [count_regressions()], [add_regression_parameters()]
#'
#' @export
extract_parameters <- function(fits, input) {

  sapply(
    names(fits),
    function(y) cbind.data.frame(
      b  = coefficients(fits[[y]]) %>%
        rprint(1),
      CI = confint.lm(fits[[y]]) %>%
        as.data.frame() %>%
        mutate_all(rprint, 1) %>%
        mutate(CI = paste0("[", `2.5 %`, ", ", `97.5 %`, "]")) %>%
        select(CI),
      t  = summary(fits[[y]])$coefficients[, "t value"]  %>%
        rprint(2),
      p  = summary(fits[[y]])$coefficients[, "Pr(>|t|)"] %>%
        sapply(function(i) pprint(p = i, text = FALSE))
    ) %>%
      rownames_to_column("x") %>%
      pivot_wider(names_from = x, values_from = -x, names_vary = "slowest")
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("y") %>%
    mutate(
      type = sapply(
        seq_len(nrow(.)),
        function(i) case_when(
          y[i] == "activities"                                                                         ~ "all",
          y[i] == "mental"   | y[i] %in% unique(with(input, subset(map, type == "mental"  )$category)) ~ "mental",
          y[i] == "physical" | y[i] %in% unique(with(input, subset(map, type == "physical")$category)) ~ "physical"
        )
      ),
      category = sapply(
        seq_len(nrow(.)),
        function(i) ifelse(
          test = y[i] %in% c("activities", "mental", "physical"),
          yes  = "all",
          no   = unique(input$map[input$map$category == y[i], "category"])
        )
      ),
      .before = 1
    ) %>%
    select(-y)
}


#' Join regression parameters onto a descriptive activity-count table
#'
#' @description
#' Left-joins the formatted regression parameters produced by
#' [extract_parameters()] onto a descriptive table produced by
#' [add_comparisons()], matching rows on the `type` and `category` columns.
#' List-columns produced by [dplyr::left_join()] are unlisted in-place.
#'
#' @param table A data frame as returned by [add_comparisons()].
#' @param fits A named list of `lm` fit objects as returned by
#'   [count_regressions()].
#' @param input A named list as returned by [get_data()], passed through to
#'   [extract_parameters()].
#'
#' @return A data frame equivalent to `table` extended with one group of
#'   four columns per regression predictor (see [extract_parameters()] for
#'   column naming conventions).
#'
#' @importFrom dplyr left_join mutate_all
#'
#' @seealso [extract_parameters()], [add_comparisons()]
#'
#' @export
add_regression_parameters <- function(table, fits, input) {
  left_join(
    table,
    extract_parameters(fits = fits, input = input),
    by = c("type", "category")
  ) %>%
    mutate_all(~ unlist(.x, use.names = FALSE))
}


#' Build binary participation contingency tables by SuperAging status
#'
#' @description
#' For each activity category, dichotomises participation (any vs. none) and
#' counts the number of SuperAgers and non-SuperAgers who participated. The
#' result also carries the activity-type label (`"mental"` / `"physical"`) and
#' per-row percentage strings for each group.
#'
#' @param input A named list as returned by [get_data()], used for
#'   `map$category` and `map$type`.
#' @param data A wide-format data frame as returned by [pivot_data()],
#'   containing the binary-count category columns and the `SA` column.
#'
#' @return A grouped data frame (by `Type`) with columns:
#'   \describe{
#'     \item{`Category`}{Activity category name.}
#'     \item{`Type`}{Activity type (`"mental"` or `"physical"`).}
#'     \item{`SA`, `nonSA`}{Raw participation counts.}
#'     \item{`SA_perc`, `nonSA_perc`}{Count and percentage strings, e.g.
#'       `"12 (60%)"`.}
#'   }
#'
#' @importFrom dplyr select mutate across group_by summarise_all ungroup
#' @importFrom tibble column_to_rownames rownames_to_column
#'
#' @seealso [chi_squares()]
#'
#' @export
contingency_tables <- function(input, data) {

  with(input, {
    data %>%
      select(SA, all_of(unique(map$category))) %>%
      mutate(across(
        .cols = all_of(unique(map$category)),
        .fns  = ~ if_else(.x > 0, TRUE, FALSE)
      )) %>%
      group_by(SA) %>%
      summarise_all(list(N = sum, perc = mean)) %>%
      column_to_rownames("SA") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("Activity") %>%
      mutate(
        Category = sub("_[^_]*$", "", Activity),
        Quant = sub("^.*_", "", Activity)
      ) %>%
      select(-Activity) %>%
      pivot_wider(
        names_from = Quant,
        values_from = c(SA, nonSA)
      ) %>%
      mutate(
        Type = unlist(
          x = sapply(
            X   = seq_len(n()),
            FUN = \(i) with(map, unique(type[category == Category[i]]))
          ),
          use.names = FALSE
        ),
        SA_perc = paste0(SA_N, " (", rprint(100 * SA_perc, 0), "%)"),
        nonSA_perc = paste0(nonSA_N, " (", rprint(100 * nonSA_perc, 0), "%)")
      ) %>%
      select(
        Category,
        Type,
        SA = SA_N,
        nonSA = nonSA_N,
        SA_perc,
        nonSA_perc
      )
  })
}


#' Run chi-square tests and compute Cramér's *V* for activity participation
#'
#' @description
#' For each activity type (`"physical"` and `"mental"`), subsets the
#' contingency table produced by [contingency_tables()] to rows with non-zero
#' totals, runs [stats::chisq.test()], and computes Cramér's *V* via
#' [rcompanion::cramerV()].
#'
#' @param tabs A data frame as returned by [contingency_tables()], containing
#'   at minimum the columns `Type`, `SA`, and `nonSA`.
#'
#' @return A data frame with one row per activity type and columns:
#'   \describe{
#'     \item{`type`}{Activity type (`"physical"` or `"mental"`).}
#'     \item{`chisq`}{Formatted chi-square statistic (3 d.p.).}
#'     \item{`df`}{Formatted degrees of freedom (0 d.p.).}
#'     \item{`p`}{Formatted *p*-value.}
#'     \item{`cramer_v`}{Formatted Cramér's *V* (3 d.p.).}
#'   }
#'
#' @importFrom stats chisq.test
#' @importFrom rcompanion cramerV
#'
#' @seealso [contingency_tables()]
#'
#' @export
chi_squares <- function(tabs) {

  sapply(
    c("physical", "mental"),
    function(x) {

      contab <- tabs %>%
        mutate(sum = SA + nonSA) %>%
        filter(Type == x & sum > 0) %>%
        select(ends_with("SA")) %>%
        as.matrix()

      chsq <- chisq.test(contab)
      v    <- cramerV(contab)

      return(data.frame(
        type     = x,
        chisq    = rprint(chsq$statistic, 3),
        df       = rprint(chsq$parameter, 0),
        p        = pprint(chsq$p.value, 3, text = FALSE),
        cramer_v = rprint(v, 3)
      ))
    }
  ) %>%
    t() %>%
    as.data.frame() %>%
    mutate_all(unlist, use.names = FALSE)
}
