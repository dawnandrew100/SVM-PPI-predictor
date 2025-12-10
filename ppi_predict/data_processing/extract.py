import pandas as pd
import re


def main():
    raw_df = pd.read_csv("../../data/BIOGRID-ALL-5.0.250.mitab.txt", sep="\t")
    print(raw_df["Interaction Types"].value_counts())
    params = ["Alt IDs Interactor A", "Alt IDs Interactor B"]
    id_dict = generate_uniprot_id_dict(raw_df, params, "Interaction Types")
    df_processed = pd.DataFrame(id_dict)
    df_processed = df_processed.rename(
        columns={
            "Alt IDs Interactor A": "Uniprot IDs Interactor A",
            "Alt IDs Interactor B": "Uniprot IDs Interactor B",
        }
    )
    df_processed["is_physical_association"] = df_processed["Interaction Types"].apply(
        lambda x: True if x == "physical association" else False
    )
    print(df_processed.head())
    print(df_processed.describe())
    print(df_processed["Interaction Types"].value_counts())
    print(df_processed["is_physical_association"].value_counts())

    df_processed.to_csv("../../data/processed.csv", index=False)


def extract_uniprot_ids(alt_id: str) -> str | None:
    """Extracts Uniprot Ids from raw string in dataframe"""

    r"""
    Positive lookbehind finds all characters equivalent to
    [a-zA-z0-9_] that have "swiss-prot:" immediately precending them.
    This only gets uniprot ids since "|" is not matched by \w
    """
    uniprot_id = re.search(r"(?<=swiss-prot:)\w+", alt_id)
    if uniprot_id:
        return uniprot_id.group()
    return None


def extract_interaction_type(interaction_type: str) -> str | None:
    """Extracts Interaction Types from raw string in dataframe"""

    """
    This regex uses a positive lookbehind to find where there's an instance
    of "( before the target matches. It then lazily matches all characters
    up to the occurence of either an opening parenthesis or closing parenthesis.
    Since it is lazy, it only matches up to the first occurence of a parenthesis.
    The parentheses are found using positive lookaheads and are included in a
    non-capturing group (?:)
    """
    interaction = re.search(r"(?<=\"\().+?(?:(?=\()|(?=\)))", interaction_type)
    if interaction:
        return interaction.group().lstrip('"')
    return None


def generate_uniprot_id_dict(
    dataframe: pd.series[str], fields: list[str], comparator: str
) -> dict[str, list[str]]:
    # Dictionary to hold fields (Interactor A and Interactor B) and Interaction Type
    data = {}
    for field in fields:
        data[field] = []
    data[comparator] = []

    for _, entry in dataframe.iterrows():
        ids = []  # Holds extracted uniprot ids per entry
        for field in fields:
            uniprot_id = extract_uniprot_ids(entry[field])
            if uniprot_id is not None:
                ids.append(uniprot_id)
        if len(ids) == len(fields):
            # Append interactors and interaction type is both ids present
            for i, field in enumerate(fields):
                data[field].append(ids[i])
            interaction = extract_interaction_type(entry[comparator])
            data[comparator].append(interaction)
    return data


if __name__ == "__main__":
    main()
