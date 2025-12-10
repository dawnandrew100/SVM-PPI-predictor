from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import requests
import json
import pandas as pd

"""
This script runs successfully HOWEVER it takes about 3-6 hours to fetch
all 58,475 sequences.
Only 58,367 sequences are actually saved when script is complete.
"""


def main():

    df = pd.read_csv("../../data/processed.csv")
    row_count = df.shape[0]
    print(df.head())

    fields = ["Uniprot IDs Interactor A", "Uniprot IDs Interactor B"]
    sequences = {}

    # Multiprocessing speeds up requests from ~1 per sec to ~4 per sec
    pool = Pool(processes=(cpu_count()))
    all_accns = pd.concat([df[field] for field in fields]).unique()
    print(f"Total unique UniProt accessions to fetch: {len(all_accns)}")

    results = pool.imap_unordered(fetch_sequence, all_accns)

    for accn, sequence in tqdm(
        results, total=len(all_accns), desc="Fetching Sequences"
    ):
        if sequence is not None:
            sequences[accn] = sequence

    pool.close()
    pool.join()

    with open("../../data/sequences.json", "w") as file:
        json.dump(sequences, file, indent=4)


def fetch_sequence(accn: str) -> (str, str | None):
    params = {"fields": ["accession", " sequence"]}
    headers = {"accept": "application/json"}
    base_url = f"https://rest.uniprot.org/uniprotkb/{accn}"

    response = requests.get(base_url, headers=headers, params=params)
    if not response.ok:
        response.raise_for_status()
        return (accn, None)

    data = response.json()
    try:
        return (accn, data["sequence"]["value"])
    except KeyError:
        return (accn, None)


if __name__ == "__main__":
    main()
