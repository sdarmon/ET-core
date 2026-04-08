# Extented-t-core

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


## *De novo* extented-t-cores computing

Currently, one can either directly use the binary version or recompile all the code and dependancies for optimal execution (recommanded).

### Dependencies and versions used (for optimal execution)

- **Python** version 3.11.2 (with **pip** v. 23.0.1 that will install the packages stored in the _requirements.txt_ file)
- **Cargo** version 1.75.0 (for Rust compilation)
- **gcc** version 12.2.0 (for C++ compilation)
- **FastP** version 0.23.4
- **BCALM 2** version v2.2.3, git commit cf371b6 (Using gatb-core version 1.4.2)
<!-- - **seqtk** version 1.4-r122 -->


### Code example (with dependencies built)

To execute the code, simply run the following command in the terminal:
```
bash extended-t-core.sh \
    --reads1 reads_1.fastq[.gz] \
    --reads2 reads_2.fastq[.gz] \
    -O output_dir
```
Where `--reads1` and `--reads2` are the paths to the paired-end reads, possibly `.gz`, and `-O` is the output directory.
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

### Binary version (without dependencies -> soon be replaced with a docker image)

```
bash extented-t-core_binaries.sh \
    --reads1 reads_1.fastq[.gz] \
    --reads2 reads_2.fastq[.gz] \
    -O output_dir
```

## Transposable Elements analysis

### Dependency and versions used

- **DFAM database** version 38 
- **Bowtie2** version 2.2.4
- **Samtools** version 1.9
- **TECount** version 1.0.0 (from the **TEtools** package)
- **featureCounts** version 2.1.1
- **Bedtools** version 2.31.1

### Code example (with dependencies built)

Will be released soon.

### Binary version (without dependencies -> soon be replaced with a docker image)

```
bash te_analysis.sh \
    --te-cons dfam_te_consensus.fa \
    -O output_dir
```
where `--te-cons` is the path to the TE consensus sequences in FASTA format, and `-O` is the output directory.
Some additional parameters can be specified, such as the number of threads to use (`-p`) and the k-mer size for the de Bruijn graph construction (`-k`). `-O`, `-p` and `-k` should be the same as the ones used for the de novo extended-t-core computing.


### Quick recap of the steps

1. Build the TE library from the **DFAM** database for **Bowtie2** alignment.
2. Align the reads to the TE library using **Bowtie2** and keep only the primary alignments.
3. Align every extended-t-core to the TE library using **Bowtie2** and keep only the primary alignments.
4. Annotate the extended-t-cores with the TE family they align to.
5. Align all the unitigs to the TE library using **Bowtie2** and keep only the primary alignments.
6. Annotate all the unitigs with the TE family they align to.
7. Compute the TE count for each TE family using **TECount**.
8. Compute the TE not covered by the extended-t-cores.
9. Save ROC curves for the TE prediction by the extended-t-cores.

### Output and files structure

The output directory will have the following additional files and directories:

    OUTDIR/
    ├── results/
    │   ├── TE_coverage_count_ab_filtered.txt       <-- Ground truth TEs
    │   ├── output_roc_curves.png                   <-- KEY VISUAL
    │   │
    │   ├── alignment/
    │   │   ├── READS.sam 
    │   │   ├── READS_sorted.bam (sorted and indexed, only primary)
    │   │   └── bowtie2_output_reads.txt
    │   │
    │   └── cores/
    │       ├── alignment_all_unitigs.sam/bam
    │       ├── all_unitigs_annotated.nodes
    │       ├── extended_t_cores_${i}_annotated.nodes
    │       └── alignment_${i}/
    │           ├── extended_t_cores_${i}_aligned.sam
    │           └──extended_t_cores_${i}_aligned.bam (sorted and indexed, only primary)
    └── data/
        └── count_TE_TECOUNT.txt


### Scripts to extract the TE consensus from the DFAM database

#### Quering of the DFAM database online

The DFAM database can be queried online to extract the TE consensus sequences for a given species 
by going to the following link : https://dfam.org/browse?classification=root%25253BInterspersed_Repeat%25253BTransposable_Element&clade_ancestors=true&clade_descendants=true

Then, you just need to specify the species name in the **Taxon** field, 
to get the list of its TE consensus sequences. One can download the results
as a fasta file, using the **FASTA** bouton at the bottom of the page.

With this link, the following options should be checked :
- Classification : *Interspersed_Repeat;Transposable_Element*
- Ancestors : *Checked*
- Descendants : *Checked*

#### Quering of the DFAM database with Curl API

The DFAM database can be queried online to extract the TE consensus sequences for a given species. 
This can be done using curl command. More details on the API can be found here : https://dfam.org/releases/Dfam_3.8/apidocs/



#### Local extraction of the TE consensus from the DFAM database

The DFAM database can be downloaded and queried locally to extract the TE consensus sequences for any given species.

Download the partitions of the **Dfam** database corresponding the studied
species (read thier README) into a `${LIBRARY_DIR}` : https://www.dfam.org/releases/current/families/FamDB/

Get **FamDB** from their GitHub : https://github.com/Dfam-consortium/FamDB

Install the Python3 package **h5py** or activate the venv (`venv/bin/activate`) of the
extended-t-core project.

```
  famdb.py -i ${LIBRARY_DIR}/ families \
  --include-class-in-name \
  --curated \
  --descendants \
  --ancestors \
  "${SPE_NAME}" --format fasta_name
```
Where Parameters :
- `-i` : path to the famdb installation
- `families` : command to extract the families
- `--include-class-in-name` : include the class of the TE in the name
- `--curated` : only curated families
- `--descendants` : include descendants of the specified species
- `--ancestors` : include ancestors of the specified species
- `${SPE_NAME}` : species name to specify (e.g. "Mus musculus")
- `--format fasta_name` : output format with fasta header containing the TE name and description