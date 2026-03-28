# SVM-PPI-predictor

Final Project for BINF 760

## Running Scripts

Python scripts can be run using `uv run <file name>.py` from the same directory of
the file.

Rust scripts can be run using `cargo run --bin <file name>` from any directory
that is within the `ppi_predict/` folder.

## Pipeline details

### Data Pre-processing

The following scripts are in the `data_processing` folder.

1. `extract.py` extracts the raw csv from the BioGRID file and creates a new dataframe
composed of the parsed uniprot IDs for protein A and protein B in the interaction.
2. `get_sequence.py` fetches protein sequences from the Uniprot database based
on their accession numbers and creates a JSON file withe each element represented
as `accession number: sequences`.
3. `filter_df_sequences` removes rows from the dataframe that do not have an
associated sequence from the previous step due to an error fetching the sequence.

### Feature addition

4. There are two scripts for feature addition: `src/bin/calc_svm_params.rs` and `src/bin/calc_paper_params.rs`.
The former is a simplified version of the latter. In both versions, the amino acid sequences
of all the proteins are pulled and converted into a triplet representation. The original alphabet
is further broken down into six groups roughly divided into shared characteristics (eg non-polar vs polar).
This allowed the alphabet to shrink from 20 to 6. These categories were converted to a numerical representation
of 0 - 5 in order to make later computations easier. The triplet representation combined with the
reduced alphabet allowed all sequences te be represented by a fixed vector with a length of 216.

This vector was converted into a relative frequency vector using the below equation.
$d_i=e^{\frac{fi-min(f1,...,f216)}{max(f1,...,f216) - min(f1,...,f216)}}-1.$

The full version (`calc_paper_params.rs`) goes a step further by also calculating the NeedlemanWunsch distance
for each sequence pair in order to add an additional parameter for calculation.

### Model

5. The final SVM mdoel is created using the previously calculated features. The model is divided into an 80% train
20% test split before being fit onto the X data (features) and the y data (whether a physical association
is present or not).

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
