# ET-core

This repository contains the code for **ET-core**, a tool designed for the _de novo_ computation of **extended-t-cores** from 
compacted de Bruijn graphs built using short-read RNA-seq data.

The **extended-t-cores** conceptually correspond to the dense regions
of the compacted de Bruijn graph due to the presence of inexact repeated sequences
in the transcriptome. There are defined as the maximal subgraphs of the compacted de
Bruijn graph where all the unitigs have a high **extended degree**. The **extended degree** of a unitig is 
a generalization of the degree of a graph node. It represents the number of locally distinct transcripts containing that node.

There are two ways to run **ET-core**: using the provided **Docker image** or by **cloning this repository**. 
While the source code version requires installing dependencies, it includes a small example dataset and scripts for
reproducibility. 
Please refer to the relevant section below based on your choice.

## ET-core from a docker image

First check if docker is installed on your machine.

If `docker ps` raise an error, you should add `sudo` to the following docker commands.

### Download the image
```
docker pull sdarmon/et-core:1.0
```
### Execution of TE-core using docker
```
docker run --rm \
    -v /absolute/path/to/your/reads/directory:/data \
    -v /absolute/path/to/your/output/directory:/output \
    sdarmon/te-core:1.0 \
    --reads1 /data/reads_1.fastq[.gz] \
    --reads2 /data/reads_2.fastq[.gz] \
    -O /output
```

Where `reads_1.fastq[.gz]` and `reads_2.fastq[.gz]` are the paired-end reads, possibly `.gz`.

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


## ET-core from the Git clone

First, clone the git project:

```
git clone https://github.com/sdarmon/ET-core
cd ET-core
```

### Dependencies and versions used (for optimal execution)

Check or install those following dependencies.

- **Python** version 3.11.2, with **pip** (version 23.0.1) to install :
    * numpy (version 2.3.1)
    * pysam (version 0.23.0)
- **Cargo** version 1.75.0 (for Rust compilation)
- **gcc** version 12.2.0 (for C++ compilation)
- **libomp-dev** version 1:18.0 (for C++ parallelisation)
- **BCALM 2** version 2.2.3, git commit cf371b6 (include gatb-core version 1.4.2)
- **Optional : FastP** version 0.23.4 (by default, run FastP on the reads. Use the `--no-fastp` option if you do NOT want to use FastP)
- **Optional : seqtk** version 1.4-r122 (for sampling the reads using the `--sample` option)


For each `<package>` ou `<xyz>` python3 package, you can use the following commands :
```
sudo apt install <package>
sudo apt install python3-<xyz> for `<xyz>` python3 package
```
### Code example (with dependencies built)

To execute the code, simply run the following command in the terminal:
```
bash ET-core.sh \
    --reads1 reads_1.fastq[.gz] \
    --reads2 reads_2.fastq[.gz] \
    -O output_dir
```
Where `reads_1.fastq[.gz]` and `reads_2.fastq[.gz]` are the paths to the paired-end reads, possibly `.gz`, and 
where `output_dir` is the output directory.

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
