DATA_URL = "https://economics.mit.edu/files/2853"


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
