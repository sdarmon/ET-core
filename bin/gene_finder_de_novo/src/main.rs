//The goal of this program is to find the genes in a graph of unitigs, starting from
//a set of unitigs that are known to be in the same component.
//Output: The list of the genes associated with the component and the subgraph
//of the neighborhood of the component.


#![allow(non_snake_case)]

use std::collections::HashMap;
use ordered_float::NotNan;

fn revcomp(seq: &str) -> String {
    seq.chars().rev().map(|c| {
        match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => c, // If the character is not A, T, C or G, we keep it as it is
        }
    }).collect()
}

//Define a function that concatenate two sequences with a overlap of K-1 bases but that can be F or R depending on the way of the edge between the two unitigs
fn concatenate(seq1: &str, seq2: &str, way1 : bool, way2 : bool, K: usize) -> String {
    if way1 && way2 {
        //FF
        return format!("{}{}", seq1, &seq2[K-1..]);
    } else if way1 && !way2 {
        //FR
        let revcomp_seq1 = revcomp(seq1);
        return format!("{}{}", seq2, &revcomp_seq1[K-1..]);
    } else if !way1 && way2 {
        //RF
        let revcomp_seq1 = revcomp(seq1);
        return format!("{}{}", revcomp_seq1, &seq2[K-1..]);
    } else {
        //RR
        return format!("{}{}", seq2, &seq1[K-1..]);
    }
}

fn main() {
    //First let's read the arguments
    let args: Vec<String> = std::env::args().collect();
    std::env::set_var("RUST_BACKTRACE", "1");
    use std::io::BufRead;
    use std::io::Write;
if args.len() != 8 {
        eprintln!("Usage: {} comp_prefix graph.nodes graph.edges graph.ab output_dir k nb_comps", args[0]);
        std::process::exit(1);
    //Input:
    //${WORK_DIR}/gene_finder_de_novo.exe \
    //           ${BASE_DIR}/comp \
    //           ${RESULTS_DIR}/graph/outputNodes.txt \
    //           ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
    //           ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance \
    //           ${BASE_DIR}/genes_of_comp \
    //           ${K} \
    //           ${MAXI}
    }

    let nb_comps = args[7].parse::<u32>().unwrap();
    if nb_comps == 0 {
        eprintln!("Error: The number of components must be greater than 0");
        std::process::exit(1);
    }

    let K = args[6].parse::<usize>().unwrap();

    //Read the graph.nodes file; it contains 3 columns separated by a tab:
    //1. Unitig ID (integer)
    //2. sequence (string)
    //3. Weight of the node (integer)

    //Read and store the graph.nodes file in dictionary with the unitig ID as key and the 2 others columns as value
    //Initialize the abundance to -1
    let mut nodes = std::collections::HashMap::new();
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&args[2]).unwrap();
    for result in reader.deserialize() {
        let (id, seq, weight): (u32, String, f64) = result.unwrap();
        nodes.insert(id, (seq, weight,-1.0));
    }

    //The edge file contains 3 columns separated by a tab:
    //1. Unitig ID 1 (integer)
    //2. Unitig ID 2 (integer)
    //3. Edge type (string: FF, FR, RF, RR)

    //Read and store the graph.edges file in a dictionary with the unitig ID as key and the list of the tuples (unitigs,edge type) it is connected to as value
    //edge_type is two letters (FF, FR, RF, RR) that indicates the orientation of the edge between the two unitigs
    let mut edges = std::collections::HashMap::new();
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&args[3]).unwrap();
    for result in reader.deserialize() {
        let (id1, id2, edge_type): (u32, u32, String) = result.unwrap();
        if !edges.contains_key(&id1) {
            edges.insert(id1, Vec::new());
        }
        edges.get_mut(&id1).unwrap().push((id2, edge_type));
    }

    //The abundance file contains 1 column:
    //1. Abundance of the ith unitig on the ith line (integer)

    //Read the abundance file and store it in the nodes dictionary
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&args[4]).unwrap();
    let mut line_number = 0;
    for result in reader.deserialize() {
        let (abundance,): (f64,) = result.unwrap();
        if let Some(node) = nodes.get_mut(&(line_number as u32)) {
            node.2 = abundance;
        }
        line_number += 1;
    }

    //Define a vector of comps_id of size nb_comps
    let mut comp_id_vec = Vec::with_capacity(nb_comps as usize);

    //Define a set of all_comp_id
    let mut all_comp_id = std::collections::HashSet::new();


    //The comps are in file prefix_comps + i + ".txt" for i in [0, nb_comps-1]
    for i in 0..nb_comps {
        let file_name = format!("{}{}.txt", &args[1], i);
        let file = std::fs::File::open(&file_name).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut comp_id = std::collections::HashSet::new();
        for line in reader.lines() {
            let id = line.unwrap().split('\t').next().unwrap().parse::<u32>().unwrap();
            comp_id.insert(id);
            all_comp_id.insert(id);
        }
        comp_id_vec.push(comp_id);
    }

    //Define max weight and max_id empty vector
    let mut max_weight = Vec::new();
    let mut max_id = Vec::new();
    //Loop over the comp_id_vec and compute the max weight and max_id for each component
    for comp_id in &comp_id_vec {
        let mut max_w = 0.0;
        let mut max_i = 0;
        for idc in comp_id {
            let (_,w,_) = nodes.get(&idc).unwrap();
            if *w > max_w {
                max_w = *w;
                max_i = *idc;
            }
        }
        max_weight.push(max_w);
        max_id.push(max_i);
    }
    //Set of the expressed transcript secondary_transcript
    let mut exp_transcript = std::collections::HashSet::new();
    let mut secondary_transcript = std::collections::HashSet::new();

    //Set of unitigs connecting components
    let mut connecting_unitigs = std::collections::HashSet::new();

    //Define a stack for the DFS using a binary heap for efficiency
    #[derive(Eq, PartialEq)]
    // depth, id, way, old_ab, start_way, initial_comp, seq_of_path
    struct StackElem(u32, u32, bool, NotNan<f64>, bool, usize, String);

    //Implement Ord and PartialOrd for StackElem for the binary heap
    //It orders by increasing depth, then by id (i.e. the first element to be poped is the lower depth)
    //Careful : BinaryHeap is a max-heap, so we need to invert the order in the cmp function
    impl Ord for StackElem {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering { //True when self.0 is less than other.0
            // Decreasing order by the first element (depth)
            other.0.cmp(&self.0) //True when other.0 is greater than self.0
                .then_with(|| self.1.cmp(&other.1)) // For tie-breaking, compare by id
        }
    }
    impl PartialOrd for StackElem {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    let mut stack = std::collections::BinaryHeap::new();

    //Add the heaviest unitig of each component to the stack
    for (i, comp_id) in comp_id_vec.iter().enumerate() {
        if !comp_id.is_empty() {
            let id = max_id[i];
            let seq = nodes.get(&id).unwrap().0.clone();
            stack.push(StackElem(0, id, true, NotNan::new(0.0).unwrap(), true, i.try_into().unwrap(),seq.clone()));
            stack.push(StackElem(0, id, false, NotNan::new(0.0).unwrap(), false, i.try_into().unwrap(), seq));
        }
    }

    //Set of the unitigs that are already visited
    let mut visited = std::collections::HashSet::new();

    //Vector of visited nodes for each component
    let mut visited_in_comp = Vec::new();
    for _ in 0..nb_comps {
        visited_in_comp.push(std::collections::HashSet::new());
    }

    //Dictionary of the comp where the unitig have been visited : visited_loc[id] = Vec<f64>
    let mut visited_loc = std::collections::HashMap::new();

    //Define the flux for each component
    let mut flux = Vec::new();
    let mut leaves = Vec::new();
    for _ in 0..nb_comps {
        flux.push(0.0);
        leaves.push(0);
    }

    //Looping over the stack
    while !stack.is_empty() {
        //Pop the first element of the stack and its way
        let StackElem(d, id, way, o_a, start_way, i,seq) = stack.pop().unwrap();
        let mut old_ab = o_a.into_inner(); // o_a était NotNan<f64]
        let mut depth = d;
        let actual_ab: f64 = nodes.get(&id).unwrap().2;
        //If the unitig is not visited
        //Check if the unitig is in the component
        let seen_id_way = visited.contains(&(id, way));
        let seen_id_in_comp = visited_in_comp[i].contains(&id);
        let seen_id_other_way = visited.contains(&(id, !way));

        if !seen_id_way && (seen_id_in_comp || !seen_id_other_way){
            visited.insert((id, way)); //Add the unitig to the visited set
            visited_in_comp[i].insert(id); //Add the unitig to the visited_in_comp set
            visited_loc.entry(id).or_insert_with(Vec::new).push((i,way,seq.clone()));//Add the unitig to the visited_loc set

            //If the unitig is in the component
            if comp_id_vec[i].contains(&id) {
                depth = 0;
                old_ab = nodes.get(&id).unwrap().2; //Get the abundance of the unitig
            }

            //A valid neighbor is an adjacent node such that its edge type is FR/FF if way == true
            //or RF/RR if way == false
            let mut neighbors = Vec::new();

            if edges.contains_key(&id) {
                for (nb, edge_type) in edges.get(&id).unwrap() {
                    if (way && (edge_type == "FR" || edge_type == "FF")) || (!way && (edge_type == "RF" || edge_type == "RR")) {
                        let new_way: bool;
                        if edge_type == "FF" || edge_type == "RF" {
                            new_way = true;
                        } else {
                            new_way = false;
                        }
                        neighbors.push((*nb, new_way));
                    }
                }
            }
            //If the unitig is a leaf
            if neighbors.len() == 0 {
                flux[i] += actual_ab;
                if old_ab > 0.0 {
                    leaves[i]+= 1;
                    exp_transcript.insert((id, depth, old_ab.to_bits(), actual_ab.to_bits(), i));
                } else {
                    secondary_transcript.insert((id, depth, actual_ab.to_bits(),i));
                }
            }

            //Add the neighbors of the unitig to the stack
            let d : u32 = (nodes.get(&id).unwrap().0.len() - K + 1).try_into().unwrap(); //Number of new unitigs added by the neighbor
            let mut max_ab = 0.0;
            let mut best_el = 0;
            for (nb, _) in &neighbors {
                let ab_n: f64 = nodes.get(&nb).unwrap().2;
                if ab_n > max_ab {
                    max_ab = ab_n;
                    best_el = *nb;
                }
            }
            for (nb, new_way) in neighbors {
                //stack.push(StackElem(depth + d, nb, new_way,  actual_ab, start_way, i.try_into().unwrap()));
                if (old_ab > 0.0) && (nb == best_el) {
                    //If the abundance of the neighbor is the greater abundance of the current node
                    //we push it to the stack with the new way, depth + d and max_ab
                    stack.push(StackElem(
                        depth + d,
                        nb,
                        new_way,
                        NotNan::new(actual_ab).expect("actual_ab must not be NaN"),
                        start_way,
                        i.try_into().unwrap(),
                        concatenate(&seq, &nodes.get(&nb).unwrap().0, way, new_way, K)
                    ));
                } else {
                    //Otherwise we push it to the stack with the new way, depth + d and 0
                    stack.push(StackElem(
                        depth + d,
                        nb,
                        new_way,
                        NotNan::new(0.0).unwrap(),
                        start_way,
                        i.try_into().unwrap(),
                        concatenate(&seq, &nodes.get(&nb).unwrap().0, way, new_way, K)
                    ));
                }
            }
        }
        if seen_id_way {
            //Case where the node has already been visited

            //First case : it has been visited from another component
            if !seen_id_in_comp {
                //Add the node to the visited_in_comp set
                visited_in_comp[i].insert(id);
                visited_loc.entry(id).or_insert_with(Vec::new).push((i,way,seq.clone()));
                //Add the node to the connecting_unitigs set
                connecting_unitigs.insert(id);
                //flux[i] += old_ab;
            } else {
                //Second case : it has been visited from the same component
                //We have a bubble
            }
        }

        if !seen_id_in_comp && seen_id_other_way { //Id seen in the other way and from another component
            //If the node has not been visited in the current component
            visited_loc.entry(id).or_insert_with(Vec::new).push((i,way,seq.clone()));
            //Add the node to the connecting_unitigs set
            connecting_unitigs.insert(id);
            //flux[i] += old_ab;
        }
    }

    //Compute for each component the number of other components it is connected to
    //by connecting_unitigs
    let mut comp_connections = vec![std::collections::HashSet::new(); nb_comps as usize];
    for id in &connecting_unitigs {
        let comp_ids = visited_loc.get(id).unwrap();
        //Add the pairwise combinations of comps in comp_ids to the comp_connections set
        for i in 0..comp_ids.len() {
            for j in (i + 1)..comp_ids.len() {
                let comp_i = comp_ids[i].0;
                let comp_j = comp_ids[j].0;
                if comp_i != comp_j {
                    comp_connections[comp_i as usize].insert(comp_j);
                    comp_connections[comp_j as usize].insert(comp_i);
                }
            }
        }
    }


    //output the results
    //Output the expressed transcript to a file

    let output_dir = &args[5];
    let mut file = std::fs::File::create(format!("{}/expressed_transcripts.txt", output_dir)).unwrap();
    for (id, depth, old_abundance, ac_ab, i) in exp_transcript {
        //Write the id, depth and abundance in the file, separated by a tab
        let seq = nodes.get(&id).unwrap().0.clone();
        writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}", id, i, seq, depth,old_abundance, ac_ab).unwrap();
    }
    //close the file
    drop(file);

    file = std::fs::File::create(format!("{}/secondary_transcripts.txt", output_dir)).unwrap();
    for (id, depth, ac_ab, i) in secondary_transcript {
        //Write the id, depth and abundance in the file, separated by a tab
        let seq = nodes.get(&id).unwrap().0.clone();
        writeln!(file, "{}\t{}\t{}\t{}\t{}", id, i, seq, depth, ac_ab).unwrap();
    }
    //close the file
    drop(file);

    //Output the connecting unitigs to a file
    file = std::fs::File::create(format!("{}/connecting_unitigs.txt", output_dir)).unwrap();
    for id in &connecting_unitigs {
        //Write the id in the file
        let seq = nodes.get(&id).unwrap().0.clone();
        let comp_ids = visited_loc.get(id).unwrap();
        let comp_ids_str = comp_ids.iter().map(|(comp_id,_,_)| comp_id.to_string()).collect::<Vec<String>>().join(",");
        let comp_path_str = comp_ids.iter().map(|(_,way, seq)| format!("{}:{}", if *way { "F" } else { "R" }, seq)).collect::<Vec<String>>().join(","); //String of the form "F:ATGC,R:GCAT"
        writeln!(file, "{}\t{}\t{}\t{}", id, seq, comp_ids_str,comp_path_str).unwrap();
    }

    //close the file
    drop(file);

    //transcript_summary_comps.txt
    file = std::fs::File::create(format!("{}/transcript_summary_comps.txt", output_dir)).unwrap();
    //Write the number of expressed transcripts for each component
    for (i, comp_id) in comp_id_vec.iter().enumerate() {
        let mut max_ab = 0.0;
        for id in comp_id {
            let ab = nodes.get(&id).unwrap().2;
            if ab > max_ab {
                max_ab = ab;
            }
        }
        writeln!(file, "Component {}: {:.2} flux \t {} max ab \t {} leaves \t {} connected comp", i, flux[i], max_ab, leaves[i], comp_connections[i].len()).unwrap();
    }
    //close the file
    drop(file);

}
