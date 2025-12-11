use serde::{Deserialize, Serialize};

use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::PathBuf;

use std::collections::BTreeMap;
use std::f64::consts::E;

// REMINDER - BTreeMap is the order by key ORD version of HashMap

// -- JSON DATA--
// BTreeMap for initial JSON data
type RawSequenceMap = BTreeMap<String, String>;

#[derive(Debug, Deserialize, Serialize)]
struct Metrics {
    sequence: String,
    rel_freq: FreqMetrics,
}

#[derive(Debug, Deserialize, Serialize)]
struct FreqMetrics {
    rel_freq: Vec<f64>,
}

// BTreeMap for final JSON data
type SequenceData = BTreeMap<String, Metrics>;


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

    let output_path = "../data/sequence_metrics.json";
    println!("Writing sequence metrics to {}...", output_path);
    write_json_file(output_path, &freq_seq_map)?;
    println!("Finished writing JSON file!");

    Ok(())
}

fn open_file(path_to_file: &str) -> Result<BufReader<File>, io::Error> {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let file_path = PathBuf::from(manifest_dir).join(path_to_file);
    let file = File::open(&file_path)?;
    let reader = BufReader::new(file);

    Ok(reader)
}

fn read_json_file(reader: BufReader<File>) -> Result<RawSequenceMap, Box<dyn Error>> {
    let raw_sequences: RawSequenceMap = serde_json::from_reader(reader)?;
    Ok(raw_sequences)
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

fn write_json_file(
    path_to_file: &str,
    data: &SequenceData,
) -> Result<(), Box<dyn Error>> {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let file_path = PathBuf::from(manifest_dir).join(path_to_file);

    if let Some(parent) = file_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let file = File::create(&file_path)?;
    serde_json::to_writer_pretty(file, data)?;

    Ok(())
}
