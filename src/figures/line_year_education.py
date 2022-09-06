import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import typer


def aggregate_data(data, cohort_limits):
    """Filters data to cohorts and aggregates data by birth quarter.

    Args:
        data (pandas.DataFrame): data to aggregate
        cohort (list[str]): cohorts to include

    Returns:
        pandas.DataFrame
    """
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


def save_line_plot(fig, path, width, height, dpi):
    """Saves a line plot of education versus quarter of birth.

    Args:
        fig (matplotlib.figure.Figure): figure to save
        path (str): path to save plot
    """
    fig.set_size_inches(width, height)
    fig.savefig(path, dpi=dpi)


def main(input_data: str, output_path: str,
         cohort: tuple[int, int] = (20, 49),
         width: int = 6, height: int = 4, dpi: int = 300):
    """Main function."""
    data = pd.read_csv(input_data)
    aggregated_data = aggregate_data(data, cohort)
    fig, ax = create_line_plot(aggregated_data)
    save_line_plot(fig, output_path, width, height, dpi)


if __name__ == "__main__":
    typer.run(main)
