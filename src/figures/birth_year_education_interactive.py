import sys

import pandas as pd
import altair as alt
from birth_year_education import aggregate_data


def create_line_plot(data):
    """Creates a line plot of education versus quarter of birth.

    Args:
        data (pandas.DataFrame): data to plot

    Returns:
        alt.Chart: line plot
    """

    brush = alt.selection_interval()

    chart_line = (
        alt.Chart(data)
        .mark_line(color="black")
        .encode(
            x=alt.X("year_quarter_of_birth:Q", title="Year of Birth"),
            y=alt.Y("education:Q", title="Education", scale=alt.Scale(zero=False)),
        )
        .properties(width=400, height=400)
    )

    chart_points = (
        alt.Chart(data)
        .mark_point(filled=True, size=80)
        .encode(
            x=alt.X("year_quarter_of_birth:Q", title="Year of Birth"),
            y=alt.Y("education:Q", title="Education"),
            color=alt.condition(
                brush,
                alt.Color("quarter_of_birth:O", title="Quarter of birth"),
                alt.value("lightgray"),
            ),  # type: ignore
        )
        .properties(width=400, height=400)
        .add_selection(brush)
    )

    return (chart_line + chart_points), brush


def create_avg_bars(data):
    """Creates a bar plot of education versus quarter of birth.

    Args:
        data (pandas.DataFrame): data to plot

    Returns:
        alt.Chart: bar plot
    """
    chart = (
        alt.Chart(data)
        .mark_bar()
        .encode(
            x=alt.X("quarter_of_birth:O", title="Quarter of birth"),
            y=alt.Y(
                "mean(education):Q", title="Education", scale=alt.Scale(zero=False)
            ),
            color=alt.Color("quarter_of_birth:O", title="Quarter of birth"),
        )
        .properties(width=200, height=400)
    )

    return chart


def create_combined_plot(data):
    """Creates a combined plot of education versus year/quarter of birth.

    Args:
        data (pandas.DataFrame): data to plot

    Returns:
        alt.Chart: combined plot
    """
    chart_line, brush = create_line_plot(data)
    chart_bars = create_avg_bars(data).transform_filter(brush)

    return chart_line | chart_bars


def main(input_data: str, output_path: str, cohort: tuple[int, int] = (20, 49)):
    data = pd.read_csv(input_data)
    aggregated_data = aggregate_data(data, cohort)
    chart = create_combined_plot(aggregated_data)
    chart.save(output_path)


if __name__ == "__main__":
    with open(snakemake.log[0], "w") as logfile:  # type: ignore # noqa: F821
        sys.stderr = sys.stdout = logfile

        main(
            input_data=snakemake.input["file"],  # type: ignore # noqa: F821
            output_path=snakemake.output["json"],  # type: ignore # noqa: F821
            cohort=snakemake.params["cohort"],  # type: ignore # noqa: F821
        )
