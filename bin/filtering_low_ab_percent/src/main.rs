#![allow(non_snake_case)]
use std::io::{Write, BufWriter};
use std::collections::HashSet;
use std::fs::File;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    std::env::set_var("RUST_BACKTRACE", "1");

    if args.len() != 6 {
        eprintln!("Usage: {} graph.nodes graph.edges graph.ab output_prefix percentage", args[0]);
        std::process::exit(1);
    }

    let percentage = args[5].parse::<f32>().expect("Percentage must be a float");
    if percentage <= 0.0 {
        eprintln!("Error: The percentage must be greater than 0");
        std::process::exit(1);
    }

    // 1. Read abundance file into a Vec (Fastest lookup)
    let mut nodes_abundance = Vec::new();
    let mut ab_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(&args[3])
        .unwrap();

    for result in ab_reader.deserialize() {
        let (abundance,): (f32,) = result.unwrap();
        nodes_abundance.push(abundance);
    }

    // 2. Prepare output files with BufWriters for speed
    let output_dir = &args[4];

    let file_edges = File::create(format!("{}_C{:.2}.edges", output_dir, percentage / 100.0)).unwrap();
    let mut writer_edges = BufWriter::new(file_edges);

    let file_removed = File::create(format!("{}_removed_C{:.2}.edges", output_dir, percentage / 100.0)).unwrap();
    let mut writer_removed = BufWriter::new(file_removed);

    let mut disconnected_nodes: HashSet<u32> = HashSet::new();

    // 3. Process edges in a single pass
    let mut edge_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(&args[2])
        .unwrap();

    for result in edge_reader.deserialize() {
        let (id1, id2, edge_type): (u32, u32, String) = result.unwrap();

        // Get abundance using usize casting
        let ab1 = nodes_abundance[id1 as usize];
        let ab2 = nodes_abundance[id2 as usize];

        // Comparison logic
        let mut is_removed = false;
        if (ab1 * 100.0) < (ab2 * percentage) {
            disconnected_nodes.insert(id1);
            is_removed = true;
        } else if (ab2 * 100.0) < (ab1 * percentage) {
            disconnected_nodes.insert(id2);
            is_removed = true;
        }

        if is_removed {
            writeln!(writer_removed, "{}\t{}\t{}", id1, id2, edge_type).unwrap();
        } else {
            writeln!(writer_edges, "{}\t{}\t{}", id1, id2, edge_type).unwrap();
        }
    }

    // Explicitly flush buffers
    writer_edges.flush().unwrap();
    writer_removed.flush().unwrap();

    // 4. Output the remaining nodes (the "Component")
    let file_nodes = File::create(format!("{}_remaining.nodes", output_dir)).unwrap();
    let mut writer_nodes = BufWriter::new(file_nodes);

    let mut node_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(&args[1])
        .unwrap();

    for result in node_reader.deserialize() {
        let (id, seq): (u32, String) = result.unwrap();

        // Only write if the node was NOT disconnected
        if !disconnected_nodes.contains(&id) {
            let ab = nodes_abundance[id as usize];
            writeln!(writer_nodes, "{}\t{}\t{}", id, seq, ab).unwrap();
        }
    }
    writer_nodes.flush().unwrap();

}