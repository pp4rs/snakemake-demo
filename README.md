# Snakemake demonstration
Reproduces some tables and graphs from Angrist &amp; Krueger (1991).

## How to compile
The project is set up so that snakemake handles the installation of the required dependencies into a local virtual environment. The following external dependencies have to be manually installed:

 * A TeX distribution with `pdflatex` and `latexmk` available on the path
 * `unrar` and `wget` available on the path
 * The `snakemake` workflow management system

Software under the first two bullet points can be installed using your preferred method. It is recommended to install snakemake in its own separate conda virtual environment (e.g. `conda create -c conda-forge -c bioconda -n snakemake snakemake`).

The steps to build the project are described in its snakemake file. If snakemake is installed it can be compiled from scratch by running the snakemake command in its root directory:

```bash
    cd /path/to/project-for-pp4rs
    conda activate snakemake
    snakemake --cores all --use-conda
```
assuming that snakemake is available in the conda environment names snakemake. `--cores all` sets the number of parallel jobs equal to the number of your logical cpu cores. If you wish to run `N` jobs in parallel, replace it with `--cores N`.
