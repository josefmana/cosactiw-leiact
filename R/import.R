#' Import and pre-process raw RETROS study data
#'
#' @description
#' Reads several sheets from an Excel workbook, harmonises participant IDs,
#' applies factor coding, and returns a named list of pre-processed data
#' frames ready for analysis. Only participants present in the SuperAging
#' classification file (`sfile`) are retained, keeping this dataset aligned
#' with the COBRA analysis sample.
#'
#' The following components are constructed:
#' \describe{
#'   \item{`dem`}{Demographic variables (age, education, physical-activity
#'     category, FAQ, MMSE).}
#'   \item{`cog`}{Cognitive test scores and SuperAging (`SA`) classification.}
#'   \item{`act`}{Long-format retrospective leisure-activity records with
#'     intensity, seasonality, and activity-category columns.}
#'   \item{`cobra`}{COBRA-A total, frequency, and time summary scores.}
#'   \item{`vls`}{VLS-AQL item-level responses in long format.}
#'   \item{`map`}{Activity mapping table linking individual activities to
#'     their type (`"mental"` / `"physical"`) and category.}
#' }
#'
#' @param mfile A character string giving the path to the activity-mapping
#'   Excel file (`.xlsx`). The mapping is read from the sheet named
#'   `"Sheet 1"` and must contain at least the columns `type`, `category`,
#'   and `activity`.
#' @param dfile A character string giving the path to the main data Excel
#'   workbook (`.xlsx`). The following sheets are expected:
#'   \describe{
#'     \item{`"1ID_DATE_AGE_SCREENING_Demogr"`}{Demographic data.}
#'     \item{`"3COGNITIVE_TESTS"`}{Cognitive test scores.}
#'     \item{`"5RETROS-LEISURE_ACTIVITIES"`}{Retrospective activity records.}
#'     \item{`"COSACTIW-pro-JAMOVI"`}{COBRA-A and VLS-AQL scores.}
#'   }
#' @param sfile A character string giving the path to an `.rds` file
#'   containing a data frame with at least the columns `id` (character,
#'   participant IDs) and `sa` (logical, SuperAging classification).
#'
#' @return A named list with elements `dem`, `cog`, `act`, `cobra`, `vls`,
#'   and `map` (see Description for details).
#'
#' @importFrom openxlsx read.xlsx
#' @importFrom dplyr select filter mutate rename left_join
#' @importFrom tidyr pivot_longer
#'
#' @seealso [pivot_data()]
#'
#' @export
get_data <- function(mfile, dfile, sfile) {

  # included women
  included <- readRDS(sfile)$id

  # activities mapping
  map <-
    read.xlsx(mfile, sheet = "Sheet 1") %>%
    select(type, category, activity) %>%
    add_row(type = "physical", category = "flexibility_health_exercise", activity = "chi_kung") %>%
    arrange(type, category, activity)

  # demography
  dem <-
    read.xlsx(dfile, sheet = "1ID_DATE_AGE_SCREENING_Demogr") %>%
    rename(
      "Age_years"            = "Age",
      "Education_years"      = "Number_of_years_of_study",
      "Education_level"      = "Highest_education_level",
      "Subjective_age_years" = "Subjective_age_number",
      "Valid"                = "Was_the_examination_valid?"
    ) %>%
    mutate(
      ID = gsub(" ", "", ID),  # no spaces in ID allowed
      PA = factor(
        PA,
        levels = 1:6,
        labels = c(
          "I'm an athletic person",
          "I enjoy movement/exercise",
          "I exercised at least 3 times a week",
          "I don't avoid movement/exercise",
          "I'm not an athletic person",
          "I had to stop doing sports (Injury)"
        )
      ),
      Education_level = factor(
        Education_level,
        levels  = 1:4,
        labels  = c(
          "Primary school",
          "Vocational school (OU)",
          "Secondary school",
          "College/university"
        ),
        ordered = TRUE
      )
    ) %>%
    filter(ID %in% included) %>%
    select(ID, Valid, Age_years, Subjective_age_years, Education_years, Education_level, PA, FAQ, MMSE)

  # cognition
  cog <-
    read.xlsx(dfile, sheet = "3COGNITIVE_TESTS") %>%
    left_join(
      readRDS(sfile) %>%
        mutate(ID = gsub(" ", "", id)) %>%
        select(ID, sa)
    ) %>%
    filter(ID %in% included) %>%
    mutate(
      ID = gsub(" ", "", ID),
      SA = factor(
        if_else(sa == TRUE, 1, 2),
        levels = 1:2,
        labels = c("SA", "nonSA")
      ),
      WHO_PA = factor(
        `WHO-PA`,
        levels = 0:1,
        labels = c("nonPA", "PA")
      )
    ) %>%
    select(ID, SA, WHO_PA, `RAVLT_1-5`, RAVLT_delayed_recall, TMT_A_time, TMT_B_time, Spon_sem, VF_animals)

  # activities
  act <-
    read.xlsx(dfile, sheet = "5RETROS-LEISURE_ACTIVITIES") %>%
    mutate(
      ID            = gsub(".0", "", gsub(" ", "", ID), fixed = TRUE),
      Activity      = if_else(Activity == "theathre", "theatre", Activity),
      Activity_type = case_when(
        Activity == "reading"        ~ "mental",
        Activity == "theatre"        ~ "mental",
        Activity == "self_education" ~ "mental",
        .default = Activity_type
      ),
      Seasonal  = if_else(Intensity >= 10, TRUE, FALSE),
      Intensity = factor(
        Intensity %% 10,  # modulo 10 to pool seasonal activities
        levels = c(0:5),
        labels = c(
          "do not know",
          "occasionally or not at all",
          "several times a month",
          "once a week",
          "several times a week",
          "every day or nearly every day"
        )
      ),
      Category = unlist(
        sapply(1:nrow(.), function(i) map[map$activity == Activity[i], "category"]),
        use.names = FALSE
      )
    ) %>%
    filter(ID %in% included)

  # COBRA-A
  cobra <-
    read.xlsx(dfile, sheet = "COSACTIW-pro-JAMOVI") %>%
    select(1, `Total-mental-activities`, `total-MA-frequency`, `Total-MA-time`) %>%
    filter(ID %in% included) %>%
    rename(
      "cobra_a_total"     = "Total-mental-activities",
      "cobra_a_frequency" = "total-MA-frequency",
      "cobra_a_time"      = "Total-MA-time"
    )

  # VLS-AQL
  vls <-
    read.xlsx(dfile, sheet = "COSACTIW-pro-JAMOVI") %>%
    select(1, starts_with("VLS_I")) %>%
    filter(ID %in% included) %>%
    pivot_longer(
      cols      = -ID,
      names_to  = "item",
      values_to = "response"
    )

  # return the data
  return(list(
    dem   = dem,
    cog   = cog,
    act   = act,
    cobra = cobra,
    vls   = vls,
    map   = map
  ))
}


#' Pivot activity data to wide format and merge with demographics and cognition
#'
#' @description
#' Takes the long-format activity records from [get_data()] and pivots them so
#' that each row is one participant and each activity, activity type, and
#' activity category is a separate column containing the binary participation
#' indicator (1 = participated, 0 = did not participate). Marginal sums across
#' all activities, activity types, and activity categories are appended, and
#' the result is joined with the cognitive and demographic data frames.
#'
#' @param input A named list as returned by [get_data()], containing at
#'   minimum the elements `act` (long-format activity data), `map` (activity
#'   mapping table), `cog` (cognitive data), and `dem` (demographic data).
#'
#' @return A wide-format data frame with one row per participant. Columns
#'   include:
#'   \describe{
#'     \item{`ID`}{Participant identifier.}
#'     \item{One column per activity name}{Binary (0/1) participation
#'       indicator.}
#'     \item{`activities`}{Total number of activities participated in.}
#'     \item{One column per activity type}{Count of activities in that type.}
#'     \item{One column per activity category}{Count of activities in that
#'       category.}
#'     \item{Columns from `cog` and `dem`}{Merged by `"ID"`.}
#'   }
#'
#' @importFrom dplyr mutate across rowSums left_join everything
#' @importFrom tibble rownames_to_column
#'
#' @seealso [get_data()]
#'
#' @export
pivot_data <- function(input) {

  with(input, sapply(
    map$activity,
    function(i) sapply(
      unique(act$ID),
      function(j) if_else(
        condition = i %in% act[act$ID == j, "Activity"],
        true  = 1,
        false = 0
      )
    )
  ) %>%
    as.data.frame() %>%
    mutate(
      activities = rowSums(across(everything())),
      !!!setNames(rep(NA, length(unique(map$type))),     unique(map$type)),
      !!!setNames(rep(NA, length(unique(map$category))), unique(map$category)),
      across(unique(map$type),     ~ rowSums(across(all_of(map[map$type     == cur_column(), "activity"])))),
      across(unique(map$category), ~ rowSums(across(all_of(map[map$category == cur_column(), "activity"]))))
    ) %>%
    rownames_to_column("ID") %>%
    left_join(cog, by = "ID") %>%
    left_join(dem, by = "ID"))
}
