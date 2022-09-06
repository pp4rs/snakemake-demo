library(optparse)
library(lfe)
library(tibble)
library(readr)
library(dplyr)
library(purrr)
library(jsonlite)
library(stringr)


load_and_filter_data <- function(data_path, cohort_limits) {
  df <- read_csv(data_path) %>%
    filter(year_of_birth >= cohort_limits[1] & year_of_birth <= cohort_limits[2])
  return(df)
}


estimate_model <- function(df, form) {
    model <- felm(form, df)
    return(model)
}


load_formula <- function(formula_path) {
    formula_string <- readLines(formula_path)
    return(as.formula(formula_string))
}


main <- function() {

    option_list <- list(
        make_option(c("--data_path"), type = "character"),
        make_option(c("--formula_path"), type = "character"),
        make_option(c("--model_out_path"), type = "character"),
        make_option(c("--cohort"), type = "character")
    )
    opt <- parse_args(OptionParser(option_list = option_list))

    limits <- as.integer(str_split(opt$cohort, ",")[[1]])
    df <- load_and_filter_data(opt$data_path, limits)
    form <- load_formula(opt$formula_path)

    model <- estimate_model(df, form)
    saveRDS(model, opt$model_out_path)

}


main()
