import pandas as pd
import json


def main():
    df = pd.read_csv("../../data/processed.csv")
    with open("../../data/sequences.json", "r") as file:
        sequence_dict = json.load(file)

    first_item = next(iter(sequence_dict.items()))
    print(first_item)
    print(df.describe())

    filtered_df = df[
        (df["Uniprot IDs Interactor A"].isin(sequence_dict.keys()))
        & (df["Uniprot IDs Interactor B"].isin(sequence_dict.keys()))
    ]

    print(filtered_df.describe())
    filtered_df.to_csv("../../data/filtered.csv", index=False)


if __name__ == "__main__":
    main()
