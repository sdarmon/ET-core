# Extented-t-core

Intro TBD

## *De novo* extented-t-cores computing

Currently, one can either directly use the binary version or recompile all the code and dependancies for optimal execution (recommanded).

### Dependencies and versions used (for optimal execution)

- **Python** version 3.11.2 (with **pip** v. 23.0.1 that will install the packages stored in the _requirements.txt_ file)
- **Cargo** version 1.75.0 (for Rust compilation)
- **gcc** version 12.2.0 (for C++ compilation)
- **FastP** version 0.23.4
- **BCALM 2**,version v2.2.3, git commit cf371b6 (Using gatb-core version 1.4.2)


### Code example (with dependencies built)

To execute the code, simply run the following command in the terminal:
```
bash extended-t-core.sh \
    --reads1 reads_1.fastq \
    --reads2 reads_2.fastq \
    -O output_dir
```
Where `--reads1` and `--reads2` are the paths to the paired-end reads, and `-O` is the output directory.
Some additional parameters can be specified, such as the number of threads to use (`-p`) and the k-mer size for the de Bruijn graph construction (`-k`), or the distance of the extended degree (`-d`).


### Quick recap of the steps

1. Build the binaries (**Rust** and **C++** codes) and create a **Python3** venvironment.
2. Run **FastP** to trim the reads, detect and remove adapters, and filter out low-quality reads. Logs are saved in the `output_dir/fastp_log.html` file.
3. Do an homopolymer compression of the reads to limit A/T stretches.
4. Build the de Bruijn graph using **BCALM2**.
5. Extract the unitigs from the de Bruijn graph and filter the sequences errors edges.
6. Compute the extended degree of the unitigs.
7. Compute the threshold for the extended degree corresponding to the top 1% of the unitigs.
8. Compute the extended-t-cores of the compacted de Bruijn graph.
9. Compute a representative sequence for each extended-t-core.
10. Classify the extended-t-cores depending on their repeat content (microsatellite, A/T stretches, others).
11. Compute the induced subgraph of the pairwise connections between the extended-t-cores.
12. (TO BE DONE) Assembly of potential full-length Transposable Elements using the extended-t-cores as seeds and the pairwise connections between them.

### Output and files structure

The output directory will contain the following files and directories:

    OUTDIR/
    ├── fastp_log.html
    │
    ├── data/
    │   ├── R1/2.fastp
    │   ├── hc1/2.fastq
    │   └── graph/
    │       ├── graph.nodes
    │       ├── graph.edges
    │       ├── graph.abundance
    │       └── unitigs_extended_degree.nodes
    │
    ├── results/
    │   ├── stretchAT_cores_list.txt
    │   ├── microsatellite_cores_list.txt
    │   ├── other_cores_list.txt
    │   ├── cores/
    │   │   └── extended_t_cores_${i}.nodes
    │   └── induced_cores_subgraph/
    │       ├── transcript_summary_comps.txt
    │       └── connecting_unitigs.txt
    │
    WORK_DIR/
    ├── script.sh
    ├── requirements.txt
    ├── venv/
    └── bin/

### Binary version (without dependencies)

```
bash extented-t-core_binaries.sh \
    --reads1 reads_1.fastq \
    --reads2 reads_2.fastq \
    -O output_dir
```

## Transposable Elements analysis

### Dependency and versions used

- **DFAM database** version XX (downloaded on XX)
- **Bowtie2** version XX
- **Samtools** version XX
- **TECount** version XX (from the **TEtools** package)
- **featureCounts** version XX
- **Bedtools** version XX

### Code example (with dependencies built)

TBD

### Binary version (without dependencies)

```
bash te_analysis.sh \
    --te-cons dfam_te_consensus.fa \
    -O output_dir
```
where `--te-cons` is the path to the TE consensus sequences in FASTA format, and `-O` is the output directory.
Some additional parameters can be specified, such as the number of threads to use (-p) and the k-mer size for the de Bruijn graph construction (-k). -O, -p and -k should be the same as the ones used for the de novo extended-t-core computing.


### Quick recap of the steps

TBD

### Output and files structure

TBD

### Script to extract the TE consensus from the DFAM database

Download and pre-requis TBD


```
  ${FAMDB_BIN} -i ${LIBRARY_DIR}/famdb/ families \
  --include-class-in-name \
  --curated \
  --descendants \
  --ancestors \
  "${SPE_NAME}" --format fasta_name \
   | sed 's/#/\t/g' \
   | sed 's/ @/\t/g' \
   | sed 's/^>/>dfam_/g' > ${DFAM_FA}
```
