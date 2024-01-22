import sys

import pandas as pd
import numpy as np


def read_data(path):
    """Reads data from path and returns a pandas DataFrame.
    Also renames columns to be more descriptive.

    Args:
        path (str): path to data

    Returns:
        pandas.DataFrame
    """

    colnames = [f"v{i+1}" for i in range(27)]
    data = pd.read_csv(
        "data/raw/QOB.txt", delim_whitespace=True, header=None, names=colnames
    )
    print(f"Read {len(data)} rows")

    rename_dict = {
        "v1": "age",
        "v2": "ageq",
        "v4": "education",
        "v5": "enocent",
        "v6": "esocent",
        "v9": "log_weekly_wage",
        "v10": "married",
        "v11": "midatl",
        "v12": "mt",
        "v13": "neweng",
        "v16": "census",
        "v18": "quarter_of_birth",
        "v19": "race",
        "v20": "smsa",
        "v21": "soatl",
        "v24": "wnocent",
        "v25": "wsocent",
        "v27": "year_of_birth",
    }
    data = data.rename(columns=rename_dict)

    return data


def prepare_data(data):
    """Prepares data for analysis.

    Args:
        data (pandas.DataFrame): data to prepare

    Returns:
        pandas.DataFrame
    """

    data = data.copy()

    data["cohort"] = np.where(
        (40 <= data["year_of_birth"]) & (data["year_of_birth"] <= 49),
        "40-49",
        np.where(
            (30 <= data["year_of_birth"]) & (data["year_of_birth"] <= 39),
            "30-39",
            "20-29",
        ),
    )
    data.loc[data["census"] == 80, "ageq"] = (
        data.loc[data["census"] == 80, "ageq"] - 1900
    )
    data["ageq_squared"] = data["ageq"] ** 2
    data["year_of_birth_within_decade"] = [
        int(str(year)[-1]) for year in data["year_of_birth"]
    ]

    return data


def main():
    """Main function."""

    with open(snakemake.log[0], "w") as logfile:  # type: ignore # noqa: F821
        sys.stderr = sys.stdout = logfile

        data = read_data(snakemake.input["file"])  # type: ignore # noqa: F821
        data = prepare_data(data)
        data.to_csv(snakemake.output["file"], index=False)  # type: ignore # noqa: F821
        print(f"Exported data to {snakemake.output['file']}")  # type: ignore # noqa: F821


if __name__ == "__main__":
    main()
