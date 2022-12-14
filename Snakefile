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
        file = "data/clean/census_data.csv",
        script = "src/figures/birth_year_education_interactive.py"
    output:
        json = "out/figures/interactive_graph.json"
    log: "logs/figures/interactive_graph.log"
    shell:
        "python {input.script} {input.file} {output.json} 2> {log}"


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
        script = "src/estimation/estimate_model.R",
        table = expand(
            "out/tables/regression_table_{cohort[0]}_{cohort[1]}.tex",
            cohort = config["cohorts"]
        )


rule regression_table:
    conda: "envs/estimation.yaml"
    input:
        script = "src/tables/regression_table.R",
        models = expand(
            "out/models/{model_name}_{model_type}_{from_}_{to}.rds",
            model_name = MODELS,
            model_type = config["model_types"],
            allow_missing=True
        )
    output:
        table = "out/tables/regression_table_{from_}_{to}.tex"
    log: "logs/tables/regression_table_{from_}_{to}.log"
    shell:
        "Rscript {input.script} --models '{input.models}' --output_path {output.table} 2> {log}"


rule all_models:
    input:
        script = "src/estimation/estimate_model.R",
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
        script = "src/estimation/estimate_model.R",
        model_spec = "src/model_specs/{model_name}.json"
    output:
        file = "out/models/{model_name}_{model_type}_{from_}_{to}.rds"
    log: "logs/models/{model_name}_{model_type}_{from_}_{to}.log"
    shell:
        "Rscript {input.script} --data_path {input.file} \
                                --formula_path {input.model_spec} \
                                --model_out_path {output.file} \
                                --cohort {wildcards.from_},{wildcards.to} \
                                --model_type {wildcards.model_type} \
                                2> {log}"


rule create_figure_schooling_diff:
    conda: "envs/figures.yaml"
    input:
        file = "data/clean/census_data.csv",
        script = "src/figures/birth_year_education.py"
    output:
        file = "out/figures/barplot_schooling_diff.png"
    log: "logs/figures/barplot_schooling_diff.log"
    shell:
        "python {input.script} barplot {input.file} {output.file} \
                --cohorts 30-39 --cohorts 40-49 \
                --width 8 --height 10 --dpi 300 2> {log}"


rule all_figures_birth_educ:
    input:
        file = expand("out/figures/line_year_education_{cohort[0]}_{cohort[1]}.png",
                      cohort = config["cohorts"])


rule create_figure_birth_educ:
    conda: "envs/figures.yaml"
    input:
        file = "data/clean/census_data.csv",
        script = "src/figures/birth_year_education.py"
    output:
        file = "out/figures/line_year_education_{from_}_{to}.png"
    log: "logs/figures/line_year_education_{from_}_{to}.log"
    shell:
        "python {input.script} lineplot {input.file} {output.file} \
                --cohort {wildcards.from_} {wildcards.to} \
                --width 8 --height 6 --dpi 300 2> {log}"


rule clean_data:
    conda: "envs/data-prep.yaml"
    input:
        file = "data/raw/QOB.txt",
        script = "src/data/prepare_data.py"
    output:
        file = "data/clean/census_data.csv"
    log: "logs/data/clean_data.log"
    shell:
        "python {input.script} {input.file} {output.file} 2> {log}"


rule download_data:
    input:
        script = "src/data/download_data.sh"
    output:
        file = "data/raw/QOB.txt"
    params:
        url = config["data_url"]
    log: "logs/data/download_data.log"
    shell:
        "bash {input.script} {params.url} {output.file} 2> {log}"


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