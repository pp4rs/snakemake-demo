library(fixest)
library(modelsummary)
library(stringr)
library(optparse)
library(purrr)
library(tibble)
library(dplyr)


create_extra_row_tibble <- function(models) {

    birth_dummies = rep("Yes", length(models))
    region_dummies = map_chr(
        models,
        function(model) ifelse("soatl" %in% rownames(model$coef), "Yes", "No")
    )
    extra_rows_list <- list(
        birth_dummies = birth_dummies,
        region_dummies = region_dummies
    )
    extra_rows <- extra_rows_list %>%
        as_tibble() %>%
        t() %>%
        as_tibble() %>% 
        mutate(term = c("9 Year-of-birth dummies", "8 Region-of-residence dummies")) %>%
        select(term, everything())

    return(extra_rows)

}


main <- function() {

    option_list <- list(
        make_option(c("--models"), type = "character",
                    help = "comma-separated paths to models"),
        make_option(c("--output_path"), type = "character")
    )

    opt <- parse_args(OptionParser(option_list = option_list))

    paths <- str_split(opt$models, " ")[[1]]
    models <- map(paths, readRDS)
    names(models) <- map_chr(models, function(model) paste(str_to_upper(model$type), model$name))

    num_models <- length(models)
    alphabetical_order <- order(names(models))
    flip_order <- rep(c(num_models / 2, 0), times = num_models / 2) + rep(seq(num_models / 2), each = 2)

    models <- models[alphabetical_order][flip_order]

    extra_rows <- create_extra_row_tibble(models)

    coef_names <- c(
        "education" = "Years of education",
        "fit_education" = "Years of education",
        "race" = "Race (1 = black)",
        "smsa" = "SMSA (1 = center city)",
        "married" = "Married (1 = married)",
        "ageq" = "Age",
        "ageq_squared" = "Age-squared"
    )

    options(modelsummary_format_numeric_latex = "mathmode")
    modelsummary(
        models,
        output = opt$output_path,
        gof_omit = ".*",
        coef_map = coef_names,
        add_rows = extra_rows
    )

}


main()
