import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import typer


app = typer.Typer()


def aggregate_data(data, cohort_limits=None):
    """Filters data to cohorts and aggregates data by birth quarter.

    Args:
        data (pandas.DataFrame): data to aggregate
        cohort (list[str]): cohorts to include

    Returns:
        pandas.DataFrame
    """
    if cohort_limits:
        data = data[data["year_of_birth"].between(*cohort_limits)]
    data = data \
        .groupby(["ageq"]) \
        .aggregate({
            "education": "mean",
            "year_of_birth": "first",
            "quarter_of_birth": "first"
        }) \
        .reset_index()
    return data


def create_line_plot(data):
    """Creates a line plot of education versus quarter of birth.

    Args:
        data (pandas.DataFrame): data to plot
        path (str): path to save plot
    """
    fig, ax = plt.subplots()
    sns.lineplot(
        data=data,
        x="ageq",
        y="education",
        color="black",
        ax=ax
    )
    sns.scatterplot(
        data=data,
        x="ageq",
        y="education",
        hue="quarter_of_birth",
        alpha=1,
        s=80,
        ax=ax
    )

    ax.set_xlabel("Year of Birth")
    ax.set_ylabel("Years of completed education")
    ax.set_title("Years of education and season of birth")

    return fig, ax


def create_bar_plot(data, cohorts):
    """Creates a line plot of education versus quarter of birth.

    Args:
        data (pandas.DataFrame): data to plot
        path (str): path to save plot
    """
    data = data.copy()
    data["education_ma5"] = data["education"].rolling(window=5, center=True).mean()
    data["education_diff_ma5"] = data["education"] - data["education_ma5"]

    fig, ax = plt.subplots(len(cohorts), 1)
    for i, cohort in enumerate(cohorts):
        data_subset = data[data["year_of_birth"].between(*cohort)]
        sns.barplot(
            data=data_subset,
            x="year_of_birth",
            y="education_diff_ma5",
            hue="quarter_of_birth",
            ax=ax[i]
        )
        ax[i].set_ylabel("Schooling differential")
        ax[i].set_xlabel("")
        if i > 0:
            ax[i].get_legend().remove()
        else:
            ax[i].get_legend().set_title('Quarter of birth')

    ax[1].set_xlabel("Year of Birth")
    fig.suptitle("Season of birth and years of schooling")

    return fig, ax


def save_plot(fig, path, width, height, dpi):
    """Saves a line plot of education versus quarter of birth.

    Args:
        fig (matplotlib.figure.Figure): figure to save
        path (str): path to save plot
    """
    fig.set_size_inches(width, height)
    fig.savefig(path, dpi=dpi)


@app.command()
def lineplot(input_data: str, output_path: str,
             cohort: tuple[int, int] = (20, 49),
             width: int = 6, height: int = 4, dpi: int = 300):
    """Create a line plot of schooling attainment by birth quarter."""
    data = pd.read_csv(input_data)
    aggregated_data = aggregate_data(data, cohort)
    fig, ax = create_line_plot(aggregated_data)
    save_plot(fig, output_path, width, height, dpi)


@app.command()
def barplot(input_data: str, output_path: str,
            cohorts: list[str] = ["30-39", "40-49"],
            width: int = 6, height: int = 8, dpi: int = 300):
    """Create a barplot of schooling differential by birth quarter."""
    data = pd.read_csv(input_data)
    cohort_tuples = []
    for cohort in cohorts:
        cohort_tuples.append(tuple(map(int, cohort.split("-"))))
    aggregated_data = aggregate_data(data, cohort_limits=None)
    fig, ax = create_bar_plot(aggregated_data, cohort_tuples)
    save_plot(fig, output_path, width, height, dpi)


if __name__ == "__main__":
    app()
