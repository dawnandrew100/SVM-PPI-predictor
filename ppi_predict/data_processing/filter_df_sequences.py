import pandas as pd
import json


def main():
    df = pd.read_csv("../../data/processed.csv")
    with open("../../data/sequences.json", "r") as file:
        sequence_dict = json.load(file)

    first_item = next(iter(sequence_dict.items()))
    print(first_item)
    print(df.describe())
    print("\nClassification split before filtering")
    print(f"{df["has_physical_association"].value_counts(normalize=True)}\n")

    filtered_df = df[
        (df["Uniprot IDs Interactor A"].isin(sequence_dict.keys()))
        & (df["Uniprot IDs Interactor B"].isin(sequence_dict.keys()))
    ]

    print(filtered_df.describe())
    print("\nClassification split after filtering")
    print(filtered_df["has_physical_association"].value_counts(normalize=True))
    filtered_df.to_csv("../../data/filtered.csv", index=False)


if __name__ == "__main__":
    main()
