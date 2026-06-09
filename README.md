# Extended-t-core

This repository contains the code to *de novo* compute the **extended-t-cores** of a 
compacted de Bruijn graph built from short reads RNA-seq.

The **extended-t-cores** conceptually correspond to the dense regions
of the compacted de Bruijn graph due to the presence of inexact repeated sequences
in the transcriptome. There are defined as the maximal subgraphs of the compacted de
Bruijn graph where all the unitigs have a high **extended degree**, a generalization
of the degree of a graph node that counts the number of neighbors at distance at
most `d` nucleotides (by default, `d=10`).



From these extended-t-cores, we can identify inexact repeats of the transcriptome,
and we can propose _de novo_ Transposable Elements (TEs) candidates.

A second script is available to compare the extended-t-cores to a TE consensus library,
such as the one from the **DFAM** database (see last section for TEs extraction).


## *De novo* Extended-t-cores computing

<!-- Currently, one can either directly use the binary version or recompile all the code and dependancies for optimal execution (recommanded).-->

### Dependencies and versions used (for optimal execution)

- **Python** version 3.11.2 (with **pip** v. 23.0.1 that will install the packages stored in the _requirements.txt_ file)
- **Python3 venv** (for the Python dependencies)
- **Cargo** version 1.75.0 (for Rust compilation)
- **gcc** version 12.2.0 (for C++ compilation)
- **FastP** version 0.23.4
- **BCALM 2** version v2.2.3, git commit cf371b6 (Using gatb-core version 1.4.2)
- **Optional : seqtk** version 1.4-r122 (for sampling the reads using the `--sample` option)


### Code example (with dependencies built)

To execute the code, simply run the following command in the terminal:
```
bash extended-t-core.sh \
    --reads1 reads_1.fastq[.gz] \
    --reads2 reads_2.fastq[.gz] \
    -O output_dir
```
Where `--reads1` and `--reads2` are the paths to the paired-end reads, possibly `.gz`, and `-O` is the output directory.

Some additional parameters can be specified:
- `-p`: number of threads to use (default: 8)
- `-k`: k-mer size to use for the DGB construction (default: 41)
- `-d`: extended degree distance to use for the weighting of the nodes (default: 10)
- `-h`: hamming distance to use for the weighting of the nodes (default: 2)
- `-t`: threshold to use for the agglomeration of the nodes (default: 'precise'; options : 'sensitive' | 'precise' | t where t is a integer greater than 1)
- `-a`: Abundance minimal for keeping the k-mers (default: 2)
- `--max-memory`: max memory to use (in MBytes, default: 14000)
- `--no-fastp` : do not run fastp on the reads (not recommended if the reads are not curated)
- `--sample` : generation a sample given a sample size of fraction (default: no sampling; options : n (number of reads) | f (sample fraction, between 0 and 1))


### Quick recap of the steps

1. Build the binaries (**Rust** and **C++** codes) and create a **Python3** venvironment.
2. Run **FastP** to trim the reads, detect and remove adapters, and filter out low-quality reads. Logs are saved in the `output_dir/fastp_log.html` file.
3. Do a homopolymer compression of the reads to limit A/T stretches.
4. Build the de Bruijn graph using **BCALM2**.
5. Extract the unitigs from the de Bruijn graph and filter the sequences errors edges.
6. Compute the extended degree of the unitigs.
7. Compute the threshold for the extended degree.
8. Compute the extended-t-cores of the compacted de Bruijn graph.
9. Compute a representative sequence for each extended-t-core.
10. Classify the extended-t-cores depending on their repeat content (microsatellite, A/T stretches, others).
11. Compute the induced subgraph of the pairwise connections between the extended-t-cores, including the DNA paths between them.

### Output and files structure

The central file is `extended_t_cores_summary.tsv`. It contains the summary of every the extended-t-cores with the 
following columns:

*The first four headers characterize the cores. The other four headers provide additional important metrics for the cores.*

| Column Header | Description | Format/Units |
| :--- | :--- | :--- |
| **`Id`** | Unique identifier for the extended $t$-core, sorted by decreasing `Max_degree`. | Integer |
| **`Representative`** | Most abundant sequence of the extended $t$-core. | DNA string |
| **`Repeat_type`** | Repeat classification based on the DNA sequences of the extended $t$-core: Microsatellite (% covered), A/T or C/G stretch $(X)^{\geq 5}$. Other cases are labeled as Potential TE. | String |
| **`TE_Score`** | Confidence score for Transposable Element classification. | `0` (Low) – `3` (High) |
| **`Max_degree`** | Highest extended degree for a node within the core. | Integer |
| **`Max_abundance`** | Highest number of reads mapped to a node of the extended $t$-core. | Float |
| **`Core_connectivity`** | Number of paths connecting other extended $t$-cores to that core. | Integer |
| **`Primary_neighbour`** | ID of the neighbour having the highest number of distinct paths connecting both cores, and percentage of such connecting paths. | ID:Percentage (%) |



The output directory will contain the following files and directories:

    OUTDIR/
    ├── fastp_log.html
    │
    ├── data/
    │   ├── R1/2.fastp.gz (reads after fastp)
    │   ├── hc1/2.fastq.gz (reads after homopolymer compression)
    │   └── graph/
    │       ├── graph.nodes (id_unitig `\t` seq_unitig)
    │       ├── graph.edges (id_unitig1 `\t` id_unitig2 `\t` way)
    │       ├── graph.abundance (abundance_unitig)
    │       └── outputNodes.txt (id_unitig `\t` seq_unitig `\t` extended_degree_unitig)
    │
    ├── results/
    │   ├── extended_t_cores_summary.tsv 
    │   ├── cores/
    │   │   └── core${i}.nodes (unitigs of the ith extended-t-core)
    │   └── induced_cores_subgraph/
    │       └── connecting_paths.txt (id_core1 `\t` id_core2 `\t` path_between_cores)
    │
    WORK_DIR/
    ├── extended-t-core.sh
    ├── requirements.txt
    ├── venv/
    └── bin/
