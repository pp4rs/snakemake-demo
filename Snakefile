from os.path import splitext, dirname
from src.utils.makeutils import find_input_files


configfile: "snake_config.yaml"
MODELS = glob_wildcards("src/model_specs/{model_name}.json").model_name


rule all:
    input:
        paper = "out/paper/paper.pdf",
        presentation = "out/presentation/presentation.html"


rule presentation:
    conda: "envs/quarto.yaml"
    input:
        qmd = "src/presentation/presentation.qmd",
        figure = "out/figures/interactive_graph.json"
    output:
        html = "out/presentation/presentation.html"
    params:
        output_dir = lambda wildcards, output: dirname(output.html)
    log: "logs/presentation/presentation.log"
    shell:
        "quarto render {input.qmd} -P figure_path:{input.figure} 2> {log} && \
         rm -rf {params.output_dir}/* && \
         mv -f src/presentation/presentation.html src/presentation/presentation_files {params.output_dir}"


rule interactive_graph:
    conda: "envs/figures.yaml"
    input:
        file = "data/clean/census_data.csv"
    output:
        json = "out/figures/interactive_graph.json"
    params:
        cohort = (20, 49)
    log: "logs/figures/interactive_graph.log"
    script:
        "src/figures/birth_year_education_interactive.py"


rule paper:
    input:
        included = find_input_files("src/paper/paper.tex"),
        tex = "src/paper/paper.tex"
    output:
        pdf = "out/paper/paper.pdf"
    params:
        output_name = lambda wildcards, output: splitext(output.pdf)[0]
    log: "logs/paper/paper.pdf"
    shell:
        "latexmk -pdf -interaction=nonstopmode -jobname={params.output_name} {input.tex} 2> {log}"


rule regression_tables_all:
    input:
        table = expand(
            "out/tables/regression_table_{cohort[0]}_{cohort[1]}.tex",
            cohort = config["cohorts"]
        )


rule regression_table:
    conda: "envs/estimation.yaml"
    input:
        models = expand(
            "out/models/{model_name}_{model_type}_{from_}_{to}.rds",
            model_name = MODELS,
            model_type = config["model_types"],
            allow_missing=True
        )
    output:
        table = "out/tables/regression_table_{from_}_{to}.tex"
    log: "logs/tables/regression_table_{from_}_{to}.log"
    script:
        "src/tables/regression_table.R"


rule all_models:
    input:
        file = expand(
            "out/models/{model_name}_{model_type}_{cohort[0]}_{cohort[1]}.rds",
            model_name = MODELS,
            model_type = config["model_types"],
            cohort = config["cohorts"]
        )


rule estimate_model:
    conda: "envs/estimation.yaml"
    input:
        file = "data/clean/census_data.csv",
        model_spec = "src/model_specs/{model_name}.json"
    output:
        file = "out/models/{model_name}_{model_type}_{from_}_{to}.rds"
    log: "logs/models/{model_name}_{model_type}_{from_}_{to}.log"
    script:
        "src/estimation/estimate_model.R"


rule create_figure_schooling_diff:
    conda: "envs/figures.yaml"
    input:
        file = "data/clean/census_data.csv",
    output:
        file = "out/figures/barplot_schooling_diff.png"
    params:
        plot_type = "barplot",
        width = 8,
        height = 6,
        dpi = 300,
        cohorts = ["30-39", "40-49"]
    log: "logs/figures/barplot_schooling_diff.log"
    script:
        "src/figures/birth_year_education.py"


rule all_figures_birth_educ:
    input:
        file = expand("out/figures/line_year_education_{cohort[0]}_{cohort[1]}.png",
                      cohort = config["cohorts"])


rule create_figure_birth_educ:
    conda: "envs/figures.yaml"
    input:
        file = "data/clean/census_data.csv",
    output:
        file = "out/figures/line_year_education_{from_}_{to}.png"
    params:
        plot_type = "lineplot",
        width = 8,
        height = 6,
        dpi = 300
    log: "logs/figures/line_year_education_{from_}_{to}.log"
    script:
        "src/figures/birth_year_education.py"


rule clean_data:
    conda: "envs/data-prep.yaml"
    input:
        file = "data/raw/QOB.txt",
    output:
        file = "data/clean/census_data.csv"
    log: "logs/data/clean_data.log"
    script:
        "src/data/prepare_data.py"


rule download_data:
    conda: "envs/download.yaml"
    input:
        script = "src/data/download_data.sh"
    output:
        file = "data/raw/QOB.txt"
    params:
        url = config["data_url"]
    log: "logs/data/download_data.log"
    shell:
        "bash {input.script} {params.url} {output.file} 2> {log}"


rule build_graphs:
    input:
        "build_graphs/filegraph.pdf",
        "build_graphs/rulegraph.pdf",
        "build_graphs/dag.pdf"


rule filegraph:
    conda: "envs/graphviz.yaml"
    input:
        "Snakefile"
    output:
        "build_graphs/filegraph.pdf"
    shell:
        "snakemake --filegraph | dot -Tpdf > build_graphs/filegraph.pdf"


rule rulegraph:
    conda: "envs/graphviz.yaml"
    input:
        "Snakefile"
    output:
        "build_graphs/rulegraph.pdf"
    shell:
        "snakemake --rulegraph | dot -Tpdf > build_graphs/rulegraph.pdf"


rule dag:
    conda: "envs/graphviz.yaml"
    input:
        "Snakefile"
    output:
        "build_graphs/dag.pdf"
    shell:
        "snakemake --dag | dot -Tpdf > build_graphs/dag.pdf"
