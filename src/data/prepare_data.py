import pandas as pd


colnames = [f"v{i+1}" for i in range(27)]
data = pd.read_csv(
    "data/raw/QOB.txt",
    delim_whitespace=True,
    header=None,
    names=colnames
)

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
    "v27": "year_of_birth"
}
data = data.rename(columns=rename_dict)
data.to_csv("data/clean/census_data.csv", index=False)
