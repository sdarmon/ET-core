
## Transposable Elements analysis

In order to reproduce the TE analysis of the extended-t-cores, one can use the `te_analysis.sh` script, which will align
the reads, the extended-t-cores and all the unitigs to a TE consensus library (e.g. from the **DFAM** database), and 
compute the TE count for each TE family using **TECount**.

### Dependency and versions used

- **DFAM database** version 38
- **Bowtie2** version 2.2.4
- **Samtools** version 1.9
- **TECount** version 1.0.0 (from the **TEtools** package)
- **featureCounts** version 2.1.1
- **Bedtools** version 2.31.1

### Code example (with dependencies built)

```
bash te_analysis.sh \
    --te-cons dfam_te_consensus.fa \
    -O output_dir
```
where `--te-cons` is the path to the TE consensus sequences in FASTA format, and `-O` is the output directory of the 
extended-t-core script. The TE consensus sequences can be extracted from the **DFAM** database using the instructions
in the last section of this README.


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

The central file is `extended_t_cores_summary_SOLUTION.tsv`. It contains the summary of the extended-t-cores, with two 
additional columns : the TE family they align to (if any) and the TE count for this family. 


The output directory will have the following additional files and directories:

    OUTDIR/
    ├── results/
    │   ├── extended_t_cores_summary_SOLUTION.tsv 
    │   ├── TE_coverage_count_ab_filtered.txt       
    │   ├── output_roc_curves.png                   
    │   │
    │   ├── alignment/
    │   │   ├── READS.sam 
    │   │   ├── READS_sorted.bam (sorted and indexed, only primary)
    │   │   └── bowtie2_output_reads.txt
    │   │
    │   └── cores/
    │       ├── alignment_all_unitigs.bam
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