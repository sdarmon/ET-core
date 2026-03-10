//The goal of this program is to find the genes in a graph of unitigs, starting from
//a set of unitigs that are known to be in the same component.
//Output: The list of the genes associated with the component and the subgraph
//of the neighborhood of the component.


#![allow(non_snake_case)]
fn main() {
    //First let's read the arguments
    let args: Vec<String> = std::env::args().collect();
    std::env::set_var("RUST_BACKTRACE", "1");

    use std::io::Write;
    use std::collections::HashSet;
if args.len() != 6 {
        eprintln!("Usage: {} graph.nodes graph.edges graph.ab output_prefix percentage", args[0]);
        std::process::exit(1);
    }

    let percentage = args[5].parse::<f32>().unwrap();
    if percentage == (0.0 as f32) {
        eprintln!("Error: The percentage must be greater than 0");
        std::process::exit(1);
    }

    //Read the graph.nodes file; it contains 2 columns separated by a tab:
    //1. Unitig ID (integer)
    //2. sequence (string)

    //Read and store the graph.nodes file in dictionary with the unitig ID as key and the 2 others columns as value
    //Initialize the abundance to -1
    let mut nodes = std::collections::HashMap::new();
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&args[1]).unwrap();
    for result in reader.deserialize() {
        let (id, seq): (u32, String) = result.unwrap();
        nodes.insert(id, (seq,0.0 as f32));
    }

    //The edge file contains 3 columns separated by a tab:
    //1. Unitig ID 1 (integer)
    //2. Unitig ID 2 (integer)
    //3. Edge type (string: FF, FR, RF, RR)

    //Function to map the string (string: FF, FR, RF, RR) to 0, 1 , 2, 3
    fn edge_type_to_int(edge_type: &str) -> u8 {
        match edge_type {
            "FF" => 0,
            "FR" => 1,
            "RF" => 2,
            "RR" => 3,
            _ => panic!("Invalid edge type: {}", edge_type),
        }
    }
    //And a reverse vector
    let edge_type_to_str = vec!["FF", "FR", "RF", "RR"];

    //Read and store the graph.edges file in a dictionary with the unitig ID as key and the list of the tuples (unitigs,edge type) it is connected to as value
    //edge_type is two letters (FF, FR, RF, RR) that indicates the orientation of the edge between the two unitigs
    let mut edges = std::collections::HashMap::new();
    let mut full_edges_set = std::collections::HashSet::new();
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&args[2]).unwrap();
    for result in reader.deserialize() {
        let (id1, id2, edge_type): (u32, u32, String) = result.unwrap();
        if !edges.contains_key(&id1) {
            edges.insert(id1, Vec::new());
        }
        edges.get_mut(&id1).unwrap().push((id2, edge_type_to_int(&edge_type)));
        full_edges_set.insert((id1,id2,edge_type_to_int(&edge_type)));
    }

    //The abundance file contains 1 column:
    //1. Abundance of the ith unitig on the ith line (integer)

    //Read the abundance file and store it in the nodes dictionary
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&args[3]).unwrap();
    let mut line_number = 0;
    for result in reader.deserialize() {
        let (abundance,): (f32,) = result.unwrap();
        if let Some(node) = nodes.get_mut(&(line_number as u32)) {
            node.1 = abundance;
        }
        line_number += 1;
    }

    let mut disconnected_nodes: HashSet<u32> = HashSet::new();
    let mut removed_edges: HashSet<(u32,u32,u8)> = HashSet::new();
    //loop over the nodes and find its neighbors in the edges dictionary
    //If one neighbor has an abundance lower than percentage*node.1, add the neighbor
    // to the disconnected_nodes set
    for (id, (_seq, ab)) in &nodes {
        for (neighbor_id, edge_type) in edges.get(id).unwrap_or(&Vec::new()) {
            let neighbor_ab = nodes.get(neighbor_id).unwrap().1;
            if (neighbor_ab * (100.0 as f32)) < (*ab  * percentage) {
                disconnected_nodes.insert(*neighbor_id);
                //Also remove the edge between the two nodes
                removed_edges.insert((*id,*neighbor_id,*edge_type));
            } else if (*ab * (100.0 as f32)) < (neighbor_ab  * percentage) { //Also check the other way around
                //Since the edges are bidirectional, we need to check both directions for the edges
                //but we need to remove only one node.
                removed_edges.insert((*id,*neighbor_id,*edge_type));
            }
        }
    }

    let remain_edges: HashSet<(u32,u32,u8)> = full_edges_set.difference(&removed_edges).cloned().collect();

    //Now we can output the new graph with the remaining edges and nodes
    let output_dir = &args[4];
    let mut file = std::fs::File::create(format!("{}_C{:.2}.edges", output_dir, percentage / (100.0 as f32))).unwrap();
    for (id1, id2, edge_type) in &remain_edges {
        writeln!(file, "{}\t{}\t{}", id1, id2, edge_type_to_str[*edge_type as usize]).unwrap();
    }
    //close the file
    drop(file);

    //Removed edges
    file = std::fs::File::create(format!("{}_removed_C{:.2}.edges", output_dir, percentage as f32 / 100.0)).unwrap();
    for (id1, id2, edge_type) in &removed_edges {
        writeln!(file, "{}\t{}\t{}", id1, id2, edge_type_to_str[*edge_type as usize]).unwrap();
    }
    //close the file
    drop(file);


    file = std::fs::File::create(format!("{}_disconnected.nodes", output_dir)).unwrap();
    for id in &disconnected_nodes {
        let seq = nodes.get(&id).unwrap().0.clone();
        let ab = nodes.get(&id).unwrap().1;
        writeln!(file, "{}\t{}\t{}", id, seq, ab).unwrap();
    }
    //close the file
    drop(file);
}
