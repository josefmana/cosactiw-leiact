# ── Formatting helpers ──────────────────────────────────────────────────────────

#' Create a directory if it does not already exist
#'
#' @description
#' A lightweight wrapper around [base::dir.create()] that checks for the
#' directory's existence first, making it safe to call repeatedly without
#' error or warning.
#'
#' @param name A character string giving the path of the directory to create.
#'
#' @return Invisibly `NULL`. Called for its side effect of creating the
#'   directory on disk.
#'
#' @seealso [base::dir.create()], [base::dir.exists()]
#'
#' @examples
#' new_folder("_figures")
#' new_folder(file.path(tempdir(), "output"))
#'
#' @export
new_folder <- function(name) {
  if (!dir.exists(name)) {
    dir.create(name, recursive = TRUE)
  }
}


#' Format a numeric value as a rounded character string
#'
#' @description
#' Converts a numeric vector to a character vector with a fixed number of
#' decimal places, using [base::sprintf()] for reliable zero-padding.
#'
#' @param x A numeric vector to be formatted.
#' @param d A single non-negative integer giving the number of decimal places.
#'   Defaults to `2`.
#'
#' @return A character vector the same length as `x`, each element formatted
#'   to `d` decimal places.
#'
#' @examples
#' rprint(3.14159)        # "3.14"
#' rprint(3.14159, d = 4) # "3.1416"
#' rprint(c(1, 2.5, NA))  # "1.00" "2.50" NA
#'
#' @export
rprint <- function(x, d = 2) {
  sprintf(paste0("%.", d, "f"), round(x, d))
}


#' Format a numeric interval as a bracketed character string
#'
#' @description
#' Wraps a pair (or vector) of numeric values in square brackets to produce a
#' human-readable confidence- or credible-interval string such as
#' `"[0.12, 0.88]"`.
#'
#' @param x A numeric vector of length two (lower bound, upper bound).
#' @param d A single non-negative integer giving the number of decimal places
#'   passed to [rprint()]. Defaults to `2`.
#'
#' @return A length-one character string of the form `"[<lower>, <upper>]"`.
#'
#' @seealso [rprint()]
#'
#' @examples
#' ciprint(c(0.1234, 0.8765))   # "[0.12, 0.88]"
#' ciprint(c(-1.5, 2.3), d = 1) # "[-1.5, 2.3]"
#'
#' @export
ciprint <- function(x, d = 2) {
  paste0("[", paste(rprint(x, d), collapse = ", "), "]")
}


#' Remove the leading zero from a decimal p-value string
#'
#' @description
#' Formats a p-value to a fixed number of decimal places and strips the
#' leading zero (e.g. `"0.043"` becomes `".043"`), following APA style.
#' Values below `.001` are returned as `"< .001"`.
#'
#' @param x A single numeric value between 0 and 1.
#' @param d A single non-negative integer giving the number of decimal places.
#'   Defaults to `3`.
#'
#' @return A length-one character string.
#'
#' @seealso [pprint()], [rprint()]
#'
#' @examples
#' zerolead(0.043)  # ".043"
#' zerolead(0.0001) # "< .001"
#'
#' @export
zerolead <- function(x, d = 3) {
  ifelse(x < .001, "< .001", sub("0.", ".", rprint(x, 3), fixed = TRUE))
}


#' Format a p-value for APA reporting
#'
#' @description
#' Converts a raw p-value to an APA-style string. Values below `.001` are
#' reported as `"< .001"`. When `text = TRUE` the string is prefixed with
#' `"= "` for inline use (e.g. `"p = .043"`).
#'
#' @param p A single numeric value between 0 and 1.
#' @param dec A single non-negative integer giving the number of decimal
#'   places. Defaults to `3`.
#' @param text Logical. If `TRUE`, non-threshold values are prefixed with
#'   `"= "`. Defaults to `FALSE`.
#'
#' @return A length-one character string.
#'
#' @seealso [zerolead()]
#'
#' @examples
#' pprint(0.043)               # ".043"
#' pprint(0.043, text = TRUE)  # "= .043"
#' pprint(0.0001)              # "< .001"
#'
#' @export
pprint <- function(p, dec = 3, text = FALSE) {
  ifelse(
    test = text == TRUE,
    yes  = ifelse(p < .001, "< .001", paste0("= ", zerolead(p, dec))),
    no   = ifelse(p < .001, "< .001",              zerolead(p, dec))
  )
}


#' Format a Wilcoxon test result as an APA string
#'
#' @description
#' Takes the list returned by [stats::wilcox.test()] and produces a
#' compact APA-style string of the form `"W = <statistic>, p <p-value>"`.
#' P-values below `.001` are shown as `"< .001"`; others are stripped of
#' their leading zero and prefixed with `"= "`.
#'
#' @param w A list as returned by [stats::wilcox.test()], containing at
#'   minimum the elements `statistic` (numeric) and `p.value` (numeric).
#'
#' @return A length-one character string, e.g. `"W = 312, p = .043"` or
#'   `"W = 312, p < .001"`.
#'
#' @seealso [wtest()], [stats::wilcox.test()]
#'
#' @examples
#' \dontrun{
#' w <- wilcox.test(rnorm(20), rnorm(20))
#' w_test(w)
#' }
#'
#' @importFrom stringr str_remove
#'
#' @export
w_test <- function(w) {
  ptext <- if (w$p.value < 0.001) {
    "< .001"
  } else {
    paste0("= ", round(w$p.value, 3) |> stringr::str_remove("^0"))
  }
  paste0("W = ", w$statistic, ", p ", ptext)
}


#' Summarise a numeric vector as "central tendency ± variability"
#'
#' @description
#' Computes a chosen measure of central tendency and a chosen measure of
#' variability for a numeric vector and concatenates them into a formatted
#' string such as `"3.14 ± 0.56"`.
#'
#' @param y A numeric vector (may contain `NA`s, which are silently removed).
#' @param dec A single non-negative integer giving the number of decimal
#'   places. Defaults to `2`.
#' @param cen A character string naming the function to use for central
#'   tendency (e.g. `"mean"`, `"median"`). Defaults to `"mean"`.
#' @param var A character string naming the function to use for variability
#'   (e.g. `"sd"`, `"IQR"`). Defaults to `"sd"`.
#' @param sep A character string used to separate the two values. Defaults to
#'   `" ± "`.
#'
#' @return A length-one character string of the form
#'   `"<central> <sep> <variability>"`.
#'
#' @seealso [rprint()]
#'
#' @examples
#' cenvar(c(1, 2, 3, 4, 5))
#' cenvar(c(1, 2, NA, 4), cen = "median", var = "IQR", sep = " | ")
#'
#' @export
cenvar <- function(y, dec = 2, cen = "mean", var = "sd", sep = " \u00b1 ") {
  sapply(
    c(cen, var),
    function(fun) do.call(fun, list(y, na.rm = TRUE)) %>% rprint(d = dec)
  ) %>%
    paste(collapse = sep)
}


#' Return the upper bound of a rounded interval
#'
#' @description
#' Given a numeric value `x`, returns a formatted string whose rounded
#' representation is guaranteed to be *at least* as large as `x` — i.e. the
#' ceiling at the requested decimal precision. Useful when reporting
#' confidence-interval upper bounds so that rounding never makes the bound
#' appear tighter than it is.
#'
#' @param x A single numeric value.
#' @param dec A single non-negative integer giving the number of decimal
#'   places. Defaults to `2`.
#'
#' @return A length-one character string.
#'
#' @seealso [rprint()]
#'
#' @examples
#' ubound(0.1251)   # "0.13" (rounded up)
#' ubound(0.12)     # "0.12" (already exact)
#'
#' @export
ubound <- function(x, dec = 2) {
  ifelse(round(x, dec) > x, rprint(x, dec), rprint(x + 1 / (10 ^ dec), dec))
}


#' Capitalise the first letter of a character string
#'
#' @description
#' Converts the first character of each element of a character vector to
#' upper-case, leaving the rest of the string untouched.
#'
#' @param x A character vector.
#'
#' @return A character vector the same length as `x`.
#'
#' @examples
#' capitalise("hello world")  # "Hello world"
#' capitalise(c("foo", "bar")) # "Foo" "Bar"
#'
#' @export
capitalise <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# ── Descriptive statistics ──────────────────────────────────────────────────────

#' Compute common descriptive statistics for a numeric variable
#'
#' @description
#' Computes n (non-missing), mean, SD, median, minimum, and maximum for a
#' single numeric column, respecting any upstream [dplyr::group_by()] grouping.
#' Uses tidy evaluation so the column name can be passed unquoted.
#'
#' @param df A data frame (optionally grouped).
#' @param x \<[`tidy-eval`][dplyr::dplyr_data_masking]\> Unquoted name of the
#'   numeric column to summarise.
#'
#' @return A one-row tibble (per group) with columns `n`, `mean`, `sd`,
#'   `median`, `min`, and `max`.
#'
#' @importFrom dplyr summarise
#' @importFrom rlang enquo `!!`
#'
#' @seealso [npc()], [cenvar()]
#'
#' @examples
#' \dontrun{
#' msm(mtcars, mpg)
#'
#' mtcars |>
#'   dplyr::group_by(cyl) |>
#'   msm(mpg)
#' }
#'
#' @export
msm <- function(df, x) {
  x <- enquo(x)
  df |>
    summarise(
      n      = sum(!is.na(!!x)),
      mean   = mean(!!x,   na.rm = TRUE),
      sd     = sd(!!x,     na.rm = TRUE),
      median = median(!!x, na.rm = TRUE),
      min    = min(!!x,    na.rm = TRUE),
      max    = max(!!x,    na.rm = TRUE)
    )
}


#' Count and percentage summary for a binary (0/1) variable
#'
#' @description
#' For a binary numeric or logical column, returns the total non-missing
#' count (`n`), the count of `1`s / `TRUE`s (`n_ano`), the count of `0`s /
#' `FALSE`s (`n_ne`), and the corresponding percentages. Respects upstream
#' [dplyr::group_by()] grouping. Uses tidy evaluation so the column name can
#' be passed unquoted.
#'
#' @param df A data frame (optionally grouped).
#' @param x \<[`tidy-eval`][dplyr::dplyr_data_masking]\> Unquoted name of the
#'   binary column (0/1 integer or logical) to summarise.
#'
#' @return A one-row tibble (per group) with columns:
#'   \describe{
#'     \item{`n`}{Total non-missing observations.}
#'     \item{`n_ano`}{Count of positive (1 / TRUE) observations.}
#'     \item{`n_ne`}{Count of negative (0 / FALSE) observations.}
#'     \item{`pc_ano`}{Percentage positive (0–100).}
#'     \item{`pc_ne`}{Percentage negative (0–100).}
#'   }
#'
#' @importFrom dplyr summarise
#' @importFrom rlang enquo `!!`
#'
#' @seealso [msm()]
#'
#' @examples
#' \dontrun{
#' df <- data.frame(frag = c(0, 1, 1, 0, 1))
#' npc(df, frag)
#' }
#'
#' @export
npc <- function(df, x) {
  x <- enquo(x)
  df |>
    summarise(
      n      = sum(!is.na(!!x)),
      n_ano  = sum(!!x,           na.rm = TRUE),
      n_ne   = n - sum(!!x,       na.rm = TRUE),
      pc_ano = mean(!!x,          na.rm = TRUE) * 100,
      pc_ne  = mean(1 - !!x,      na.rm = TRUE) * 100
    )
}


# ── Inferential statistics ──────────────────────────────────────────────────────

#' Run an APA-formatted independent-samples t-test
#'
#' @description
#' Splits `df` by `group`, extracts the two groups' values of column `x`, and
#' passes them to [apa::t_test()] / [apa::t_apa()] to obtain a formatted
#' APA result string (e.g. `"t(38) = 2.10, p = .043, d = 0.66 [0.01,
#' 1.30]"`).
#'
#' @param df A data frame containing at least the columns `x` and `group`.
#' @param x A character string giving the name of the numeric outcome column.
#' @param group A character string giving the name of the grouping column.
#'   Must have exactly two distinct values. Defaults to `"frag"`.
#'
#' @return A character string with the APA-formatted t-test result, as
#'   produced by [apa::t_apa()].
#'
#' @importFrom dplyr pull
#'
#' @seealso [wtest()], [w_test()], [apa::t_apa()]
#'
#' @examples
#' \dontrun{
#' ttest(df = my_data, x = "score", group = "condition")
#' }
#'
#' @export
ttest <- function(df, x, group = "frag") {
  gvals <- sort(unique(df[[group]]))
  stopifnot(length(gvals) == 2)
  g1 <- df[df[[group]] == gvals[1], x] |> pull()
  g2 <- df[df[[group]] == gvals[2], x] |> pull()
  apa::t_apa(apa::t_test(g1, g2))
}


#' Run an independent-samples Wilcoxon rank-sum test
#'
#' @description
#' Splits `df` by `group`, extracts the two groups' values of column `x`, and
#' passes them to [stats::wilcox.test()]. Use [w_test()] to format the
#' returned object as an APA string.
#'
#' @param df A data frame containing at least the columns `x` and `group`.
#' @param x A character string giving the name of the numeric outcome column.
#' @param group A character string giving the name of the grouping column.
#'   Must have exactly two distinct values. Defaults to `"frag"`.
#'
#' @return An object of class `htest` as returned by [stats::wilcox.test()].
#'
#' @importFrom dplyr pull
#'
#' @seealso [w_test()], [ttest()], [stats::wilcox.test()]
#'
#' @examples
#' \dontrun{
#' w <- wtest(df = my_data, x = "score", group = "condition")
#' w_test(w)
#' }
#'
#' @export
wtest <- function(df, x, group = "frag") {
  gvals <- sort(unique(df[[group]]))
  stopifnot(length(gvals) == 2)
  g1 <- df[df[[group]] == gvals[1], x] |> pull()
  g2 <- df[df[[group]] == gvals[2], x] |> pull()
  wilcox.test(g1, g2)
}


# ── ggplot2 themes & tables ─────────────────────────────────────────────────────

#' APA-style ggplot2 theme
#'
#' @description
#' A [ggplot2::theme()] add-on that approximates APA 7th-edition figure
#' formatting: no background grid, solid axis lines, bold axis titles, and a
#' legend with a black border. Designed to be layered on top of
#' [ggplot2::theme_minimal()].
#'
#' @details
#' Theme adapted from
#' <https://ggplot2tutor.com/tutorials/apa_theme_rprofile>.
#'
#' @return A [ggplot2::theme()] object.
#'
#' @importFrom ggplot2 theme element_rect element_text element_line
#'   element_blank unit margin
#'
#' @seealso [ggplot2::theme_minimal()], [gt_apa()]
#'
#' @examples
#' \dontrun{
#' tibble(x = 1:10, y = 1:10) |>
#'   ggplot(aes(x = x, y = y)) +
#'   geom_point() +
#'   theme_minimal(base_size = 18) +
#'   apa_theme()
#' }
#'
#' @export
apa_theme <- function() {
  ggplot2::theme(
    plot.margin          = ggplot2::unit(c(1, 1, 1, 1), "cm"),
    plot.background      = ggplot2::element_rect(fill = "white", color = NA),
    plot.title           = ggplot2::element_text(
      size   = 22, face = "bold",
      hjust  = 0.5,
      margin = ggplot2::margin(b = 15)
    ),
    axis.line            = ggplot2::element_line(color = "black", linewidth = .5),
    axis.title           = ggplot2::element_text(size = 18, color = "black", face = "bold"),
    axis.text            = ggplot2::element_text(size = 15, color = "black"),
    axis.text.x          = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
    axis.title.y         = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
    axis.ticks           = ggplot2::element_line(linewidth = .5),
    panel.grid           = ggplot2::element_blank(),
    legend.position.inside = c(0.20, 0.8),
    legend.background    = ggplot2::element_rect(color = "black"),
    legend.text          = ggplot2::element_text(size = 15),
    legend.margin        = ggplot2::margin(t = 5, l = 5, r = 5, b = 5),
    legend.key           = ggplot2::element_rect(color = NA, fill = NA)
  )
}


#' Visualise posterior predictive distributions via \pkg{bayesplot}
#'
#' @description
#' A convenience wrapper around [bayesplot::pp_check()] that draws posterior
#' predictive samples from a fitted \pkg{brms} model and compares them to
#' observed data, optionally grouping by a categorical variable. Supports
#' density-overlay and summary-statistic plot types.
#'
#' @param fit A fitted model object accepted by [brms::posterior_predict()]
#'   (typically a `brmsfit`).
#' @param data A data frame containing both the observed outcome column and the
#'   grouping column.
#' @param y A character string giving the name of the outcome column in `data`.
#' @param x A character string giving the name of the grouping column in
#'   `data`.
#' @param labs A list of length two: `labs[[1]]` is the plot title and
#'   `labs[[2]]` is the subtitle. Either element may be `NULL`. Defaults to
#'   `list(NULL, NULL)`.
#' @param meth A character string passed as the `fun` argument to
#'   [bayesplot::pp_check()]. Defaults to `"ppc_dens_overlay_grouped"`.
#' @param draws An integer vector of draw indices to pass as the row subscript
#'   of the `yrep` matrix. Defaults to a random sample of 100 draws from
#'   1–4000.
#' @param stat A character string naming the summary statistic for
#'   `"stat_grouped"` plot types (e.g. `"mean"`, `"sd"`). Pass `NULL` for
#'   density-overlay methods. Defaults to `"mean"`.
#'
#' @return A `ggplot` object.
#'
#' @importFrom brms posterior_predict
#' @importFrom bayesplot pp_check
#' @importFrom ggplot2 labs theme element_text
#'
#' @seealso [brms::posterior_predict()], [bayesplot::pp_check()]
#'
#' @export
ppc_plot <- function(
    fit,
    data,
    y, x,
    labs  = list(NULL, NULL),
    meth  = "ppc_dens_overlay_grouped",
    draws = sample(1:4000, 100),
    stat  = "mean"
) {
  pp_check(
    object = data[, y],
    yrep   = posterior_predict(fit, newdata = data)[draws, ],
    fun    = meth,
    stat   = stat,
    group  = data[, x]
  ) +
    labs(
      title    = labs[[1]],
      subtitle = labs[[2]]
    ) +
    theme(
      legend.position = "none",
      plot.title      = element_text(hjust = .5, face = "bold"),
      plot.subtitle   = element_text(hjust = .5)
    )
}


#' Apply APA-style formatting to a \pkg{gt} table
#'
#' @description
#' Pipes a data frame through [gt::gt()] and applies a set of border, font,
#' alignment, and background options that approximate APA 7th-edition table
#' formatting (white outer borders, thick black lines above and below the
#' column headers, and a plain white table body).
#'
#' @param x A data frame or tibble to be formatted.
#' @param grp A character string giving the name of the column to use as the
#'   group label (passed to `groupname_col` in [gt::gt()]). Defaults to
#'   `NULL` (no row grouping).
#' @param nms A character string giving the name of the column to use as the
#'   stub / row-name column (passed to `rowname_col` in [gt::gt()]). Defaults
#'   to `NULL`.
#' @param title A character string (may contain HTML) for the table title.
#'   Defaults to `" "` (a single space, which suppresses the visible heading
#'   while preserving its spacing).
#'
#' @return A `gt_tbl` object ready for display or further modification.
#'
#' @importFrom gt gt tab_options px pct cols_align tab_style cell_borders
#'   cell_text cell_fill cells_body tab_header html opt_align_table_header
#'
#' @seealso [gt::gt()], [gt::tab_options()], [apa_theme()]
#'
#' @examples
#' \dontrun{
#' mtcars[1:4, 1:4] |> gt_apa(title = "Descriptive Statistics")
#' }
#'
#' @export
gt_apa <- function(x, grp = NULL, nms = NULL, title = " ") {
  x %>%
    gt(groupname_col = grp, rowname_col = nms) %>%
    tab_options(
      table.border.top.color            = "white",
      heading.title.font.size           = px(16),
      column_labels.border.top.width    = 3,
      column_labels.border.top.color    = "black",
      column_labels.border.bottom.width = 3,
      column_labels.border.bottom.color = "black",
      table_body.border.bottom.color    = "black",
      table.border.bottom.color         = "white",
      table.width                       = pct(100),
      table.background.color            = "white"
    ) %>%
    cols_align(align = "center") %>%
    tab_style(
      style = list(
        cell_borders(
          sides  = c("top", "bottom"),
          color  = "white",
          weight = px(1)
        ),
        cell_text(align = "center"),
        cell_fill(color = "white", alpha = NULL)
      ),
      locations = cells_body(
        columns = everything(),
        rows    = everything()
      )
    ) %>%
    tab_header(title = html(title)) %>%
    opt_align_table_header(align = "left")
}

