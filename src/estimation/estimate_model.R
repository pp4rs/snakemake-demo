library(optparse)
library(fixest)
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
    model <- feols(form, df)
    model$data <- NULL
    model$model <- NULL

    return(model)
}


load_formula <- function(specs, iv=FALSE) {
    dep_var <- specs$dep_var
    indep_vars <- unlist(specs$indep_vars)

    if (!is.null(specs$fixed_effects)) {
        fixed_effects_part <- unlist(specs$fixed_effects)
    } else {
        fixed_effects_part <- "0"
    }


    if (iv) {
        instrumented_vars <- unlist(specs$instrumental$instrumented_vars)
        instruments <- unlist(specs$instrumental$instruments)
        indep_vars <- setdiff(indep_vars, instrumented_vars)
        if (length(indep_vars) == 0) {
            indep_vars <- "1"
        }

        iv_part <- paste(
            paste(instrumented_vars, collapse = " | "),
            "~",
            paste(instruments, collapse = " + ")
        )
    }

    formula_str <- paste(
        dep_var, "~",
        paste(indep_vars, collapse = " + "), "|",
        fixed_effects_part
    )
    if (iv) {
        formula_str <- paste(formula_str, "|", iv_part)
    }

    return(as.formula(formula_str))
}


main <- function() {

    option_list <- list(
        make_option(c("--data_path"), type = "character"),
        make_option(c("--formula_path"), type = "character"),
        make_option(c("--model_out_path"), type = "character"),
        make_option(c("--cohort"), type = "character"),
        make_option(c("--model_type"), type = "character")
    )
    opt <- parse_args(OptionParser(option_list = option_list))

    if (opt$model_type == "iv") {
        iv <- TRUE
    } else if (opt$model_type == "ols") {
        iv <- FALSE
    } else {
        stop("model_type must be either 'iv' or 'ols'")
    }

    specs <- read_json(opt$formula_path)

    limits <- as.integer(str_split(opt$cohort, ",")[[1]])
    df <- load_and_filter_data(opt$data_path, limits)
    form <- load_formula(specs, iv=iv)

    model <- estimate_model(df, form)
    model$name <- specs$name
    model$type <- opt$model_type
    saveRDS(model, opt$model_out_path)

}


main()
