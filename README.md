# SVM-PPI-predictor

Final Project for BINF 760

## Running Scripts

Python scripts can be run using `uv run <file name>.py` from the same directory of
the file.

Rust scripts can be run using `cargo run --bin <file name>` from any directory
that is within the `ppi_predict/` folder.

## Pipeline details

### Data Pre-processing

1. `extract.py` extracts the raw csv from the BioGRID file and creates a new dataframe
composed of the parsed uniprot IDs for protein A and protein B in the interaction.
2. `get_sequence.py` fetches protein sequences from the Uniprot database based
on their accession numbers and creates a JSON file withe each element represented
as `accession number: sequences`.
3. `filter_df_sequences` removes rows from the dataframe that do not have an
associated sequence from the previous step due to an error fetching the sequence.

## Data Disclosure

The data used in this pipeline was acquired from
[BioGRID](https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-5.0.250/)
and is available for use to copy, modify, and publish under the MIT license.

The specific download used was `BIOGRID-ALL-5.0.250.mitab.zip`. Due to size
constraints, this file is not included in this repository, but can be easily
accessed from the BioGRID website.

## Citations

- Stark C, Breitkreutz BJ, Reguly T, Boucher L, Breitkreutz A, Tyers M.
BioGRID: a general repository for interaction datasets. Nucleic Acids Res. 2006 Jan 1;
34(Database issue):D535-9. doi: 10.1093/nar/gkj109. PMID: 16381927; PMCID: PMC1347471.
