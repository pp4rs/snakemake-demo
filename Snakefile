DATA_URL = "https://economics.mit.edu/files/2853"


rule create_figure_schooling_diff:
    conda: "envs/figures.yaml"
    input:
        file = "data/clean/census_data.csv",
        script = "src/figures/birth_year_education.py"
    output:
        file = "out/figures/barplot_schooling_diff.png"
    shell:
        "python {input.script} barplot {input.file} {output.file} \
                --cohorts 30-39 --cohorts 40-49 \
                --width 8 --height 10 --dpi 300"



rule all_figures_birth_educ:
    input:
        file = expand("out/figures/line_year_education_{cohort[0]}_{cohort[1]}.png",
                      cohort=[(30, 39), (40, 49)])


rule create_figure_birth_educ:
    conda: "envs/figures.yaml"
    input:
        file = "data/clean/census_data.csv",
        script = "src/figures/line_year_education.py"
    output:
        file = "out/figures/birth_year_education_{from_}_{to}.png"
    shell:
        "python {input.script} lineplot {input.file} {output.file} \
                --cohort {wildcards.from_} {wildcards.to} \
                --width 8 --height 6 --dpi 300"


rule clean_data:
    conda: "envs/data-prep.yaml"
    input:
        file = "data/raw/QOB.txt",
        script = "src/data/prepare_data.py"
    output:
        file = "data/clean/census_data.csv"
    shell:
        "python {input.script} {input.file} {output.file}"


rule download_data:
    input:
        script = "src/data/download_data.sh"
    output:
        file = "data/raw/QOB.txt"
    shell:
        "bash {input.script} {DATA_URL} {output.file}"
