#!/bin/bash
#Need some comments here


begin=`date +%s`


P="8"
K="41"
D_NT="10"
T=""
H="2"
S=""
A="2"
READS_1=""
READS_2=""
OUTDIR=""
MAX_MEM="14000"
NO_FASTP=""
VALIDE=""
# Parse options
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p)
            P="$2"
            shift 2
            ;;
        -k)
            K="$2"
            shift 2
            ;;
        -d)
            D_NT="$2"
            shift 2
            ;;
        -t)
            T="$2"
            shift 2
            ;;
        -h)
            H="$2"
            shift 2
            ;;
        -a)
            A="$2"
            shift 2
            ;;
        --reads1)
            READS_1="$2"
            shift 2
            ;;
        --reads2)
            READS_2="$2"
            shift 2
            ;;
        --sample)
            S="$2"
            shift 2
            ;;
        --max-memory)
            MAX_MEM="$2"
            shift 2
            ;;
        --no-fastp)
            NO_FASTP="True"
            SKIP_FASTP="True"
            shift 1
            ;;
        --help)
            HELP="True"
            shift 1
            ;;
        -O)
            OUTDIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            VALIDE="F"
            shift 1
            ;;
    esac
done

#If READS_1, READS_2 or OUTDIR are not set, exit
if [[ -z "$READS_1" || -z "$READS_2" || -z "$OUTDIR" || -n ${HELP} || -n ${VALIDE} ]]; then
    echo -e "Usage : $0  \n\t --reads1 <reads1.fastq[.gz]> \n\t --reads2 <reads2.fastq[.gz]> \n\t -O <output_dir> \n\t [-p <threads>] \n\t [-k <k-mer size>] \n\t [-d <extended degree distance>] \n\t [-h <hamming distance>] \n\t [-t <threshold>] \n\t [-a <abundance min>] \n\t [--max-memory <MB>]  \n\t [--no-fastp] \n\t [--sample <sample size | sample frac>]  \n\t [--help] \n"

    echo "Mandatory arguments : "
    echo -e "\t --reads1: path to the first reads file (fastq or fastq.gz)"
    echo -e "\t --reads2: path to the second pair-ended reads file (fastq or fastq.gz)"
    echo -e "\t -O: path to the output directory"
    echo -e "-> fasta and fasta.gz files are usable only if the option '--no-fastp' is toggled\n"

    echo "Optional arguments : "
    echo -e "\t -p: number of threads to use (default: 8)"
    echo -e "\t -k: k-mer size to use for the DBG construction (default: 41)"
    echo -e "\t -d: extended degree distance to use for the weighting of the nodes (default: 10)"
    echo -e "\t -h: hamming distance to use for the weighting of the nodes (default: 2)"
    echo -e "\t -t: threshold to use for the agglomeration of the nodes (default: 'precise'; options : 'sensitive' | 'precise' | t where t is a integer greater than 1)"
    echo -e "\t -a: Abundance minimal for keeping the k-mers (default: 2)\n"

    echo "Miscellaneous arguments : "
    echo -e "\t --max-memory: max memory to use (in MBytes, default: 14000)"
    echo -e "\t --no-fastp : do not run fastp on the reads (not recommended if the reads are not curated)"
    echo -e "\t --sample : generation a sample given a sample size of fraction (default: no sampling; options : n (number of reads) | f (sample fraction, between 0 and 1))"
    echo -e "\t --help: display this help message and exit \n"
    exit 1
fi

#If OUTDIR end with a /, remove it
if [[ "$OUTDIR" == */ ]]; then
    OUTDIR="${OUTDIR%/}"
fi

# Check values
echo "Threads: $P"
echo "K-mer: $K"
echo "Extended Degree distance: $D_NT"
echo "Reads1: $READS_1"
echo "Reads2: $READS_2"
echo "Output dir: $OUTDIR"

#Local variables
DATA_DIR=${OUTDIR}/data
RESULTS_DIR=${OUTDIR}/results
BASE_DIR=${RESULTS_DIR}/cores
# WORK_DIR = path to the directory of this script
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR=${WORK_DIR}/bin

#Create the general directories of that species
mkdir -p $DATA_DIR
mkdir -p $RESULTS_DIR
mkdir -p $BASE_DIR
mkdir -p $DATA_DIR/graph
mkdir -p ${RESULTS_DIR}/induced_cores_subgraph

if [[ -z "${SKIP_BUILD_RUST}" ]]; then
  ##Build the bin for the gene_finder function
  echo "Building the bin for the gene_finder function..."
  cargo build --release --manifest-path ${BIN_DIR}/gene_finder_de_novo/Cargo.toml
  cp ${BIN_DIR}/gene_finder_de_novo/target/release/gene_finder_de_novo ${BIN_DIR}/gene_finder_de_novo.exe

  ##Build the bin for the filtering_low_ab_percent function
  echo "Building the bin for the filtering_low_ab_percent function..."
  cargo build --release --manifest-path ${BIN_DIR}/filtering_low_ab_percent/Cargo.toml
  cp ${BIN_DIR}/filtering_low_ab_percent/target/release/filtering_low_ab_percent ${BIN_DIR}/filtering_low_ab_percent.exe

  ##Build the bin for the filtering_low_ab_percent function
  echo "Building the bin for the homopolymorphic compression..."
  cargo build --release --manifest-path ${BIN_DIR}/at_compressor/Cargo.toml
  cp ${BIN_DIR}/at_compressor/target/release/at_compressor ${BIN_DIR}/homopolymorphic_compression.exe

  g++ -O3 -g -fopenmp ${BIN_DIR}/graph.cpp ${BIN_DIR}/ponderation.cpp -o ${BIN_DIR}/graph.exe
  g++ -O3 -g ${BIN_DIR}/graph.cpp ${BIN_DIR}/agglo.cpp -o ${BIN_DIR}/agglo.exe

  #Removng default venv
  #python3 -m venv venv
  #source ${WORK_DIR}/venv/bin/activate
  #pip install -r requirements.txt
  #deactivate

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo -e "Build time (in seconds): $elapsed \n"
fi

  #source ${WORK_DIR}/venv/bin/activate

#Sample the reads if the sample size is specified, and write the sampled reads in ${READS_1}.sampled and ${READS_2}.sampled
if [[ -n "${S}" ]]; then

    if [[ -z "${SKIP_SAMPLING}" ]]; then
      echo "Sampling the reads"
      seqtk sample -s 42 ${READS_1} ${S} | gzip -c > ${READS_1}.${S}.sampled.gz
      seqtk sample -s 42 ${READS_2} ${S} | gzip -c > ${READS_2}.${S}.sampled.gz
    fi 
    # Update READS_1 and READS_2 to point to the sampled files
    READS_1="${READS_1}.${S}.sampled.gz"
    READS_2="${READS_2}.${S}.sampled.gz"

    #Check if the sampled files are empty, if yes, exit
    if [[ ! -s "${READS_1}" || ! -s "${READS_2}" ]]; then
        echo "Error: Sampled files are empty. Please check the sample size and the input files."
        echo "For large sample sizes, use fractional sampling (e.g., --sample 0.1 for 10% of the reads) instead of a fixed number of reads."
        exit 1
    fi
fi

#Start the benchmark time now
start=`date +%s`
begin=`date +%s`
#Skipping fastp using user's files
if [[ -n "${NO_FASTP}" ]]; then
  if [[ "${READS_1}" == *.gz ]]; then
    cp ${READS_1} ${DATA_DIR}/R1.fastp.gz
  else
    gzip -c ${READS_1} > ${DATA_DIR}/R1.fastp.gz
  fi

  if [[ "${READS_2}" == *.gz ]]; then
    cp ${READS_2} ${DATA_DIR}/R2.fastp.gz
  else
    gzip -c ${READS_2} > ${DATA_DIR}/R2.fastp.gz
  fi
fi

#Skipping fastp using fastP output (for debugging)
if [[ -z "${SKIP_FASTP}" ]]; then
  ##FastP of the reads to remove the poly(A) tails
  echo "FastP of the reads ..."

  fastp \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --trim_poly_x \
      --thread ${P} \
      --poly_x_min_len 5 \
      -z 4 \
      --html ${RESULTS_DIR}/fastp_log.html \
      --in1 ${READS_1} \
      --in2 ${READS_2} \
      --out1 ${DATA_DIR}/R1.fastp.gz \
      --out2 ${DATA_DIR}/R2.fastp.gz
fi



if [[ -z "${SKIP_HC}" ]]; then
  ##Compute the HC
  echo "HC of the reads ..."
    ${BIN_DIR}/homopolymorphic_compression.exe  ${DATA_DIR}/R1.fastp.gz ${DATA_DIR}/hc_1.fa.gz 5
    ${BIN_DIR}/homopolymorphic_compression.exe  ${DATA_DIR}/R2.fastp.gz ${DATA_DIR}/hc_2.fa.gz 5
    end=`date +%s`
    elapsed=`expr $end - $begin`
    begin=`date +%s`
    echo -e " HC (in seconds): $elapsed \n"
fi
if [[ -z "${SKIP_BCALM}" ]]; then
    ##Compute the DGB with bcalm
    echo "DBG with bcalm ..."
    #ls -1 ${DATA_DIR}/hc_1.fa ${DATA_DIR}/hc_2.fa > ${DATA_DIR}/list_reads
    echo "hc_1.fa.gz" > ${DATA_DIR}/list_reads
    echo "hc_2.fa.gz" >> ${DATA_DIR}/list_reads

  M="10"
  #If ${K} <= 10, we change M to half the value of K
  if (( $((${K} < 11)) )); then
    M=$((${K}/2))
  fi

  bcalm \
      -in ${DATA_DIR}/list_reads \
      -kmer-size ${K} \
      -nb-cores ${P} \
      -abundance-min ${A} \
      -minimizer-size ${M} \
      -max-memory ${MAX_MEM} \
      -out ${DATA_DIR}/graph/hc_1_hc_2_k${K}

  rm  ${DATA_DIR}/graph/hc_1_hc_2_k${K}.unitigs.fa.glue*

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo -e "Bcalm time (in seconds): $elapsed \n"
fi



if [[ -z "${SKIP_FILTERING}" ]]; then
awk -f ${BIN_DIR}/bcalm_unitig_to_edges.awk ${DATA_DIR}/graph/hc_1_hc_2_k${K}.unitigs.fa > ${DATA_DIR}/graph/hc_1_hc_2_k${K}.edges
awk '/^>/ {id = substr($1, 2); next } {print id "\t" $0}' ${DATA_DIR}/graph/hc_1_hc_2_k${K}.unitigs.fa > ${DATA_DIR}/graph/hc_1_hc_2_k${K}.nodes
awk '/^>/ {ab = substr($4, 6); print ab}' ${DATA_DIR}/graph/hc_1_hc_2_k${K}.unitigs.fa > ${DATA_DIR}/graph/hc_1_hc_2_k${K}.abundance
  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo -e "DBG time (in seconds): $elapsed \n"

  ##Filter the low abundance percent unitigs
  echo "Filtering the low abundance percent unitigs ..."
  ${BIN_DIR}/filtering_low_ab_percent.exe \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}.nodes \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}.edges \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}.abundance \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K} \
      5

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo -e "Filtering time (in seconds): $elapsed \n"
fi

if [[ -z "${SKIP_GEN_GRAPH}" ]]; then

  ##Compute the weighting
  echo "Weighting of the nodes..."
  ${BIN_DIR}/graph.exe \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}.nodes \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}_C0.05.edges \
      ${D_NT} \
      -k ${K}  \
      -o ${DATA_DIR}/graph/unitigs_extended_degree.nodes  \
      -h ${H}


  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo -e "Weighting time (in seconds): $elapsed \n"
fi

if [[ -z "${SKIP_THRESHOLD}" ]]; then
  ##Compute the threshold
  echo "Threshold of the nodes..."
  if [[ -z "${T}" || "${T}" == "precise" ]]; then
    # Default case or explicit sensitive
    T=$(python3 "${BIN_DIR}/plot.py" "${DATA_DIR}/graph/unitigs_extended_degree.nodes" top0001)
  elif [[ "${T}" == "sensitive" ]]; then
    T=$(python3 "${BIN_DIR}/plot.py" "${DATA_DIR}/graph/unitigs_extended_degree.nodes" top1)
  fi
 echo "T=${T}"

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo -e "Threshold time (in seconds): $elapsed \n"
fi

if [[ -z "${SKIP_AGGLO}" ]]; then
  ## Compute the connexe cores
  echo "Agglomeration of connexe cores..."
  #If the BASE_DIR is not empty, remove the files in it
  shopt -s nullglob
   for file in ${BASE_DIR}/*; do
        rm -r ${file}
      done
shopt -u nullglob

  ${BIN_DIR}/agglo.exe \
      ${DATA_DIR}/graph/unitigs_extended_degree.nodes \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}_C0.05.edges \
      -c ${T} \
      -k ${K} \
      ${BASE_DIR}


  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo -e "Agglomeration time (in seconds): $elapsed \n"

fi

if [[ -z "${SKIP_CONSENSUS}" ]]; then
  echo "Computing the representative sequences of the cores..."
    ##The number of cores to compute
    MAXI=0 #$(ls ${BASE_DIR}/core*.txt | wc -l)
    shopt -s nullglob
    for file in ${BASE_DIR}/core*.txt; do
      MAXI=$(( $MAXI + 1 ))
    done
    shopt -u nullglob
    echo "Number of cores : ${MAXI}"


    ##Compute the analysis over every core
    python3 ${BIN_DIR}/seq_consensium_of_cores.py \
        ${BASE_DIR}/core \
        ${DATA_DIR}/graph/hc_1_hc_2_k${K}_C0.05.edges \
        ${DATA_DIR}/graph/hc_1_hc_2_k${K}.abundance \
        ${RESULTS_DIR}/seq_consensium.txt \
        ${K} \
        ${MAXI}

      end=`date +%s`
      elapsed=`expr $end - $begin`
      begin=`date +%s`
      echo -e "Representative sequences time (in seconds): $elapsed \n"
fi

if [[ -z "${SKIP_INDUCED}" ]]; then
     echo "Analysis of the extended_t_cores induced subgraph..."
    ${BIN_DIR}/gene_finder_de_novo.exe \
          ${BASE_DIR}/core \
          ${DATA_DIR}/graph/unitigs_extended_degree.nodes \
          ${DATA_DIR}/graph/hc_1_hc_2_k${K}_C0.05.edges \
          ${DATA_DIR}/graph/hc_1_hc_2_k${K}.abundance \
          ${RESULTS_DIR}/induced_cores_subgraph \
          ${K} \
          ${MAXI}

     awk '{print $3}' FS='\t' ${RESULTS_DIR}/induced_cores_subgraph/connecting_unitigs.txt | sed 's/,/\t/g' > ${RESULTS_DIR}/induced_cores_subgraph/connecting_edges.txt
    python3 ${BIN_DIR}/connecting_to_edges.py ${RESULTS_DIR}/induced_cores_subgraph/connecting_edges.txt  > ${RESULTS_DIR}/induced_cores_subgraph/connected.edges



    end=`date +%s`
    elapsed=`expr $end - $begin`
    begin=`date +%s`
    echo -e "Analysis time (in seconds): $elapsed \n"
fi



if [[ -z "${SKIP_CLASSIFICATION}" ]]; then
  echo -e "De novo classification of the extended_t_cores... \n"
    python3 ${BIN_DIR}/analysis_core_de_novo.py \
        ${BASE_DIR}/core \
        ${DATA_DIR}/graph/hc_1_hc_2_k${K}.abundance \
        ${RESULTS_DIR}/seq_consensium.txt \
        ${RESULTS_DIR}/analysis_core \
        ${RESULTS_DIR}/induced_cores_subgraph/connecting_paths.txt



    printf "Core_ID\tRepresentative\tMax abundance\n" >  ${RESULTS_DIR}/microsatellite_cores_list.txt
    sort -k3,3nr ${RESULTS_DIR}/analysis_core_microsat.txt.temp >> ${RESULTS_DIR}/microsatellite_cores_list.txt
    rm ${RESULTS_DIR}/analysis_core_microsat.txt.temp

    printf "Core_ID\tRepresentative\tMax abundance\n" >  ${RESULTS_DIR}/stretchAT_cores_list.txt
    sort -k3,3nr ${RESULTS_DIR}/analysis_core_stretchAT.txt.temp >> ${RESULTS_DIR}/stretchAT_cores_list.txt
    rm ${RESULTS_DIR}/analysis_core_stretchAT.txt.temp

    printf "Core_ID\tRepresentative\tMax abundance\n" >  ${RESULTS_DIR}/stretchCG_cores_list.txt
    sort -k3,3nr ${RESULTS_DIR}/analysis_core_stretchCG.txt.temp >> ${RESULTS_DIR}/stretchCG_cores_list.txt
    rm ${RESULTS_DIR}/analysis_core_stretchCG.txt.temp

    printf "Core_ID\tRepresentative\tMax abundance\n" >  ${RESULTS_DIR}/other_cores_list.txt
    sort -k3,3nr ${RESULTS_DIR}/analysis_core_others.txt.temp >> ${RESULTS_DIR}/other_cores_list.txt
    rm ${RESULTS_DIR}/analysis_core_others.txt.temp

      end=`date +%s`
      elapsed=`expr $end - $begin`
      begin=`date +%s`
      echo -e "De novo classification time (in seconds): $elapsed \n"
fi

total_end=`date +%s`
total_elapsed=`expr $total_end - $start`
echo "Total time (in h:m:s): $(printf '%02d:%02d:%02d\n' $(($total_elapsed/3600)) $(($total_elapsed%3600/60)) $(($total_elapsed%60)))"
