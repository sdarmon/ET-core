#!/bin/bash
#Need some comments here


begin=`date +%s`

P=""
K=""
READS_1=""
READS_2=""
OUTDIR=""

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
        --reads1)
            READS_1="$2"
            shift 2
            ;;
        --reads2)
            READS_2="$2"
            shift 2
            ;;
        -O)
            OUTDIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check values
echo "Threads: $P"
echo "K-mer: $K"
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


  g++ -O3 -g graph.cpp ponderation.cpp -o graph.exe
  g++ -O3 -g ${BIN_DIR}/graph.cpp ${BIN_DIR}/agglo.cpp -o ${BIN_DIR}/agglo.exe

  python3 -m venv venv
  source venv/bin/activate
  pip install -r requirements.txt


  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Build time (in seconds): $elapsed \n"
fi


if [[ -z "${SKIP_FASTP}" ]]; then
  ##FastP of the reads to remove the poly(A) tails
  echo "FastP of the reads ..."

  FastP \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --trim_poly_x \
      --thread ${P} \
      --poly_x_min_len 5 \
      --html ${RESULTS_DIR}/fastp_log.html \
      --in1 ${READS_1} \
      --in2 ${READS_2} \
      --out1 ${DATA_DIR}/R1.fastp \
      --out2 ${DATA_DIR}/R2.fastp
fi



if [[ -z "${SKIP_HC}" ]]; then
  ##Compute the HC
  echo "HC of the reads ..."
  python3 ${BIN_DIR}/homomorphic_compression.py  ${DATA_DIR}/R1.fastp ${DATA_DIR}/hc_1.fq 5
  python3 ${BIN_DIR}/homomorphic_compression.py  ${DATA_DIR}/R2.fastp ${DATA_DIR}/hc_2.fq 5

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "FastP time and HC (in seconds): $elapsed \n"

  ##Compute the DGB with bcalm
  echo "DGB with bcalm ..."
  ls -1 ${DATA_DIR}/hc_1.fq ${DATA_DIR}/hc_2.fq > ${DATA_DIR}/list_reads
  bcalm \
      -in ${DATA_DIR}/list_reads \
      -kmer-size ${K} \
      -abundance-min 2 \
      -nb-cores ${P} \
      -out ${DATA_DIR}/graph/hc_1_hc_2_k${K}
fi


if [[ -z "${SKIP_GEN_GRAPH}" ]]; then

awk -f bcalm_unitig_to_edges.awk ${DATA_DIR}/graph/hc_1_hc_2_k${K}.unitigs.fa > ${DATA_DIR}/graph/hc_1_hc_2_k${K}.edges
awk '/^>/ {id = substr($1, 2); next } {print id "\t" $0}' ${DATA_DIR}/graph/hc_1_hc_2_k${K}.unitigs.fa > ${DATA_DIR}/graph/hc_1_hc_2_k${K}.nodes
awk '/^>/ {ab = substr($4, 6); print ab}' ${DATA_DIR}/graph/hc_1_hc_2_k${K}.unitigs.fa > ${DATA_DIR}/graph/hc_1_hc_2_k${K}.abundance
  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "DGB time (in seconds): $elapsed \n"

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
  echo "Filtering time (in seconds): $elapsed \n"


  ##Compute the weighting
  echo "Weighting of the nodes..."
  ${BIN_DIR}/graph.exe \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}.nodes \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}_C0.05.edges \
      ${D_NT} \
      -k ${K}  \
      -o ${DATA_DIR}/graph/outputNodes.txt


  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Weighting time (in seconds): $elapsed \n"


  ##Compute the threshold
  echo "Threshold of the nodes..."
  if [[ -z "${SKIP_AGGLO}" ]]; then
      T=$(python3 ${BIN_DIR}/plot.py ${DATA_DIR}/graph/outputNodes.txt top1)
      # Update the T variable in the environment.sh file
      #sed -i "s/^T=.*/T=${T}/" "${ENV}"
      echo "T=${T}"
  fi

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Threshold time (in seconds): $elapsed \n"
fi

if [[ -z "${SKIP_AGGLO}" ]]; then
  ## Compute the connexe components
  echo "Agglomeration of connexe components..."
  #rm -r ${BASE_DIR}/*
   for file in ${BASE_DIR}/*; do
        rm ${file}
      done

  ${BIN_DIR}/agglo.exe \
      ${DATA_DIR}/graph/outputNodes.txt \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}_C0.05.edges \
      -c ${T} \
      -d ${D_NT} \
      ${RESULTS_DIR} \
      -clean ${DATA_DIR}/graph/hc_1_hc_2_k${K}.abundance > ${RESULTS_DIR}/rapportAgglo.txt


  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Agglomeration time (in seconds): $elapsed \n"

fi

if [[ -z "${SKIP_CONSENSUS}" ]]; then
    ##The number of components to compute
     MAXI=0 #$(ls ${BASE_DIR}/comp*.txt | wc -l)
        for file in ${BASE_DIR}/comp*.txt; do
          MAXI=$(( $MAXI + 1 ))
        done
        echo "Number of comps : ${MAXI}"


    ##Compute the analysis over every component
    echo "Consensus sequence :" #>> ${BASE_DIR}/induced_cores_subgraph${i}/gene_summary.txt
    python3 ${BIN_DIR}/seq_consensium_of_comps.py \
        ${BASE_DIR}/comp \
        ${DATA_DIR}/graph/hc_1_hc_2_k${K}_C0.05.edges \
        ${DATA_DIR}/graph/hc_1_hc_2_k${K}.abundance \
        ${RESULTS_DIR}/seq_consensium.txt \
        ${K} \
        ${MAXI}

      end=`date +%s`
      elapsed=`expr $end - $begin`
      begin=`date +%s`
      echo "Consensus sequences time (in seconds): $elapsed \n"

    ${BIN_DIR}/gene_finder_de_novo.exe \
          ${BASE_DIR}/comp \
          ${DATA_DIR}/graph/outputNodes.txt \
          ${DATA_DIR}/graph/hc_1_hc_2_k${K}_C0.05.edges \
          ${DATA_DIR}/graph/hc_1_hc_2_k${K}.abundance \
          ${RESULTS_DIR}/induced_cores_subgraph \
          ${K} \
          ${MAXI}

     awk '{print $3}' FS='\t' ${RESULTS_DIR}/induced_cores_subgraph/connecting_unitigs.txt | sed 's/,/\t/g' > ${RESULTS_DIR}/induced_cores_subgraph/connecting_edges.txt
    python3 ${BIN_DIR}/connecting_to_edges.py ${RESULTS_DIR}/induced_cores_subgraph/connecting_edges.txt  > ${RESULTS_DIR}/induced_cores_subgraph/connected.edges


    python3 ${BIN_DIR}/analysis_comp_de_novo.py \
        ${BASE_DIR}/comp \
        ${DATA_DIR}/graph/hc_1_hc_2_k${K}.abundance \
        ${RESULTS_DIR}/seq_consensium.txt \
        ${RESULTS_DIR}/analysis_comp

      end=`date +%s`
      elapsed=`expr $end - $begin`
      begin=`date +%s`
      echo "Gene finding and analysis time (in seconds): $elapsed \n"
fi