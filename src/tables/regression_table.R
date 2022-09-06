library(lfe)
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
    names(models) <- paste("(", seq(1, length(models)), ")", sep = "")

    #extra_rows <- create_extra_row_tibble(models)

    coef_names <- c(
        "education" = "Years of education",
        "education(fit)" = "Years of education",
        "race" = "Race (1 = black)",
        "smsa" = "SMSA (1 = center city)",
        "married" = "Married (1 = married)",
        "ageq" = "Age",
        "ageq_squared" = "Age-squared"
    )

    modelsummary(
        models,
        output = opt$output_path,
        gof_omit = ".*",
        coef_map = coef_names#,
        #add_rows = extra_rows
    )

}


main()
