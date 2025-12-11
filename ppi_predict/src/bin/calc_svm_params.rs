use csv::{Reader, Writer};
use serde::{Deserialize, Serialize};

use goombay_rs::align::{AlignmentMatrices, NeedlemanWunsch};

use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::PathBuf;

use std::collections::BTreeMap;
use std::f64::consts::E;

// REMINDER - BTreeMap is the order by key ORD version of HashMap
// THIS SCRIPT TAKES ~30hrs TO RUN

// -- JSON DATA--
// BTreeMap for initial JSON data
type RawSequenceMap = BTreeMap<String, String>;

#[derive(Debug, Deserialize)]
struct Metrics {
    sequence: String,
    rel_freq: FreqMetrics,
}

#[derive(Debug, Deserialize)]
struct FreqMetrics {
    rel_freq: Vec<f64>,
}

// BTreeMap for final JSON data
type SequenceData = BTreeMap<String, Metrics>;

// -- CSV DATA --
#[derive(Debug, Deserialize)]
struct InteractionRecord {
    #[serde(rename = "Uniprot IDs Interactor A")]
    interactor_a: String,
    #[serde(rename = "Uniprot IDs Interactor B")]
    interactor_b: String,
    #[serde(rename = "Interaction Types")]
    interaction_types: String,
    has_physical_association: String,
}

#[derive(Debug, Serialize)]
struct NewInteractionRecord {
    #[serde(rename = "Uniprot IDs Interactor A")]
    interactor_a: String,
    #[serde(rename = "Interactor A Frequency Vector")]
    interactor_a_di: String,
    #[serde(rename = "Uniprot IDs Interactor B")]
    interactor_b: String,
    #[serde(rename = "Interactor B Frequency Vector")]
    interactor_b_di: String,
    #[serde(rename = "Interaction Types")]
    interaction_types: String,
    #[serde(rename = "Normalised Distance Score")]
    normalized_distance: f64,
    has_physical_association: String,
}

fn main() -> Result<(), Box<dyn Error>> {
    let reader = open_file("../data/sequences.json")?;
    let raw_sequences: RawSequenceMap = read_json_file(reader)?;

    let category_nums = ['0', '1', '2', '3', '4', '5'];
    let categories = vec![
        vec!['I', 'V', 'L', 'M'], // 0
        vec!['F', 'Y', 'W'],      // 1
        vec!['H', 'K', 'R'],      // 2
        vec!['D', 'E'],           // 3
        vec!['Q', 'N', 'T', 'P'], // 4
        vec!['A', 'C', 'G', 'S'], // 5
    ];

    let mut freq_seq_map: SequenceData = BTreeMap::new();
    let feature_map = gen_freq_map(&category_nums);
    println!("Generating frequency map for each accession...");
    for (accn, seq) in raw_sequences {
        let rel_freq = calculate_metrics(&seq, &categories, &feature_map);
        let metrics = Metrics {
            sequence: seq.clone(),
            rel_freq,
        };
        freq_seq_map.insert(accn, metrics);
    }
    println!("Finished generating sequence metrics!");

    // Output path for new CSV
    let writer_file = create_output_file("../data/svm_features.csv")?;
    let mut wtr = Writer::from_writer(writer_file);

    println!("Adding headers to file...");
    let headers: Vec<String> = vec![
        "Uniprot IDs Interactor A".to_string(),
        "Interactor A Frequency Vector".to_string(),
        "Uniprot IDs Interactor B".to_string(),
        "Interactor B Frequency Vector".to_string(),
        "Interaction Types".to_string(),
        "Normalised Distance Score".to_string(),
        "has_physical_association".to_string(),
    ];
    wtr.write_record(&headers)?;
    println!("Successfully added headers!");

    println!("Adding rows to CSV file...");
    // Add values to new columns of CSV
    let reader = open_file("../data/filtered.csv")?;
    let mut rdr = read_csv_file(reader)?;
    for result in rdr.deserialize() {
        let record: InteractionRecord = result?;

        let interactor_a_id = &record.interactor_a;
        let interactor_b_id = &record.interactor_b;

        let sequence_a: &String;
        let d_a_vector: String;
        let sequence_b: &String;
        let d_b_vector: String;

        match freq_seq_map.get(interactor_a_id) {
            Some(sequence_metrics) => {
                sequence_a = &sequence_metrics.sequence;
                d_a_vector = vec_to_string(&sequence_metrics.rel_freq.rel_freq);
            }
            None => continue,
        }
        match freq_seq_map.get(interactor_b_id) {
            Some(sequence_metrics) => {
                sequence_b = &sequence_metrics.sequence;
                d_b_vector = vec_to_string(&sequence_metrics.rel_freq.rel_freq);
            }
            None => continue,
        }

        let nw = NeedlemanWunsch::compute(sequence_a, sequence_b);
        let norm_dist = nw.normalized_distance();

        let new_record = NewInteractionRecord {
            interactor_a: interactor_a_id.clone(),
            interactor_a_di: d_a_vector,
            interactor_b: interactor_b_id.clone(),
            interactor_b_di: d_b_vector,
            interaction_types: record.interaction_types.clone(),
            normalized_distance: norm_dist,
            has_physical_association: record.has_physical_association.clone(),
        };

        wtr.serialize(new_record)?;
    }
    wtr.flush()?;
    println!("Successfully wrote all interaction records to CSV.");

    Ok(())
}

fn open_file(path_to_file: &str) -> Result<BufReader<File>, io::Error> {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let file_path = PathBuf::from(manifest_dir).join(path_to_file);
    let file = File::open(&file_path)?;
    let reader = BufReader::new(file);

    Ok(reader)
}

fn create_output_file(path_to_file: &str) -> Result<File, io::Error> {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let file_path = PathBuf::from(manifest_dir).join(path_to_file);

    if let Some(parent) = file_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    File::create(file_path)
}

fn read_json_file(reader: BufReader<File>) -> Result<RawSequenceMap, Box<dyn Error>> {
    let raw_sequences: RawSequenceMap = serde_json::from_reader(reader)?;
    Ok(raw_sequences)
}

fn read_csv_file(reader: BufReader<File>) -> Result<Reader<BufReader<File>>, Box<dyn Error>> {
    let rdr = csv::Reader::from_reader(reader);
    Ok(rdr)
}

fn calculate_metrics(
    sequence: &str,
    categories: &[Vec<char>],
    feature_map: &BTreeMap<String, usize>,
) -> FreqMetrics {
    let mut feature_map = feature_map.clone();
    let window_size = 3;
    let end = sequence.len();
    for i in 0..=(end - window_size) {
        let triplet = &sequence[i..i + window_size];
        let category_triplets = aa_to_category(triplet, categories);
        feature_map
            .entry(category_triplets.to_string())
            .and_modify(|count| *count += 1);
    }

    let f_min = feature_map
        .values()
        .cloned()
        .reduce(|a, b| if a < b { a } else { b })
        .unwrap() as f64;
    let f_max = feature_map
        .values()
        .cloned()
        .reduce(|a, b| if a > b { a } else { b })
        .unwrap() as f64;

    let denominator = f_max - f_min;
    if denominator == 0.0 {
        let rel_freq = feature_map.iter().map(|(_, _)| 0.0).collect();

        return FreqMetrics { rel_freq };
    }

    let mut freq_vec = vec![0.0; 216];
    for (i, (_, &f_i)) in feature_map.iter().enumerate() {
        let exponent = (f_i as f64 - f_min) / denominator;
        let exponential_term = E.powf(exponent);
        let d_i = exponential_term - 1.0;

        freq_vec[i] += d_i;
    }
    FreqMetrics { rel_freq: freq_vec }
}

fn gen_freq_map(terms: &[char]) -> BTreeMap<String, usize> {
    let mut all_codes: Vec<String> = Vec::new();
    for i in terms {
        for j in terms {
            for k in terms {
                let code = format!("{}{}{}", i, j, k);
                all_codes.push(code);
            }
        }
    }
    let freq_map: BTreeMap<String, usize> = all_codes.into_iter().map(|k| (k, 0_usize)).collect();

    freq_map
}

fn aa_to_category(triplet: &str, categories: &[Vec<char>]) -> String {
    let mut category_triplet = Vec::new();
    for aa in triplet.chars() {
        for (i, group) in categories.iter().enumerate() {
            if group.contains(&aa) {
                category_triplet.push(i.to_string());
                continue;
            }
        }
    }
    category_triplet.into_iter().collect()
}

fn vec_to_string(vec: &Vec<f64>) -> String {
    let temp_string = vec.iter()
    .map(|f| f.to_string())
    .collect::<Vec<String>>()
    .join(", ");

    format!("[{}]", temp_string)
}
