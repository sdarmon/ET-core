#!/bin/bash
#Need some comments here


begin=`date +%s`
start=`date +%s`


P="8"
K="41"
D_NT="10"
T=""
H="2"
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

#If READS_1, READS_2 or OUTDIR are not set, exit
if [[ -z "$READS_1" || -z "$READS_2" || -z "$OUTDIR" ]]; then
    echo "Use: $0  --reads1 <reads1.fastq[.gz]> --reads2 <reads2.fastq[.gz]> -O <output_dir> [-p <threads>] [-k <k-mer size>] [-d <extended degree distance>]"
    echo "Default values: -p 8 -k 41 -d 10 -h 2"
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

  g++ -O3 -g ${BIN_DIR}/graph.cpp ${BIN_DIR}/ponderation.cpp -o ${BIN_DIR}/graph.exe
  g++ -O3 -g ${BIN_DIR}/graph.cpp ${BIN_DIR}/agglo.cpp -o ${BIN_DIR}/agglo.exe

  python3 -m venv venv
  source ${WORK_DIR}/venv/bin/activate
  pip install -r requirements.txt
  deactivate

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Build time (in seconds): $elapsed \n"
fi


  source ${WORK_DIR}/venv/bin/activate

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
    echo "FastP time and HC (in seconds): $elapsed \n"

    ##Compute the DGB with bcalm
    echo "DGB with bcalm ..."
    #ls -1 ${DATA_DIR}/hc_1.fa ${DATA_DIR}/hc_2.fa > ${DATA_DIR}/list_reads
    echo "hc_1.fa.gz" > ${DATA_DIR}/list_reads
    echo "hc_2.fa.gz" >> ${DATA_DIR}/list_reads
  bcalm \
      -in ${DATA_DIR}/list_reads \
      -kmer-size ${K} \
      -abundance-min 2 \
      -nb-cores ${P} \
      -out ${DATA_DIR}/graph/hc_1_hc_2_k${K}
fi



if [[ -z "${SKIP_FILTERING}" ]]; then
awk -f ${BIN_DIR}/bcalm_unitig_to_edges.awk ${DATA_DIR}/graph/hc_1_hc_2_k${K}.unitigs.fa > ${DATA_DIR}/graph/hc_1_hc_2_k${K}.edges
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
fi

if [[ -z "${SKIP_GEN_GRAPH}" ]]; then

  ##Compute the weighting
  echo "Weighting of the nodes..."
  ${BIN_DIR}/graph.exe \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}.nodes \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}_C0.05.edges \
      ${D_NT} \
      -k ${K}  \
      -o ${DATA_DIR}/graph/outputNodes.txt  \
      -h ${H}


  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Weighting time (in seconds): $elapsed \n"
fi

if [[ -z "${SKIP_THRESHOLD}" ]]; then
  ##Compute the threshold
  echo "Threshold of the nodes..."
      if [[ -z "${T}" ]]; then
        T=$(python3 ${BIN_DIR}/plot.py ${DATA_DIR}/graph/outputNodes.txt top1)
      fi
      echo "T=${T}"


  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Threshold time (in seconds): $elapsed \n"
fi

if [[ -z "${SKIP_AGGLO}" ]]; then
  ## Compute the connexe components
  echo "Agglomeration of connexe components..."
  #If the BASE_DIR is not empty, remove the files in it
  shopt -s nullglob
   for file in ${BASE_DIR}/*; do
        rm -r ${file}
      done
shopt -u nullglob

  ${BIN_DIR}/agglo.exe \
      ${DATA_DIR}/graph/outputNodes.txt \
      ${DATA_DIR}/graph/hc_1_hc_2_k${K}_C0.05.edges \
      -c ${T} \
      -k ${K} \
      ${BASE_DIR}


  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Agglomeration time (in seconds): $elapsed \n"

fi

if [[ -z "${SKIP_CONSENSUS}" ]]; then
  echo "Computing the representative sequences of the components..."
    ##The number of components to compute
    MAXI=0 #$(ls ${BASE_DIR}/comp*.txt | wc -l)
    shopt -s nullglob
    for file in ${BASE_DIR}/comp*.txt; do
      MAXI=$(( $MAXI + 1 ))
    done
    shopt -u nullglob
    echo "Number of comps : ${MAXI}"


    ##Compute the analysis over every component
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
      echo "Representative sequences time (in seconds): $elapsed \n"
fi

if [[ -z "${SKIP_INDUCED}" ]]; then
     echo "Analysis of the extended_t_cores induced subgraph..."
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

    end=`date +%s`
    elapsed=`expr $end - $begin`
    begin=`date +%s`
    echo "Analysis time (in seconds): $elapsed \n"
fi



if [[ -z "${SKIP_CLASSIFICATION}" ]]; then
  echo "De novo classification of the extended_t_cores... \n"
    python3 ${BIN_DIR}/analysis_comp_de_novo.py \
        ${BASE_DIR}/comp \
        ${DATA_DIR}/graph/hc_1_hc_2_k${K}.abundance \
        ${RESULTS_DIR}/seq_consensium.txt \
        ${RESULTS_DIR}/analysis_comp


    printf "Comp_ID\tRepresentative\tMax abundance\n" >  ${RESULTS_DIR}/microsatellite_cores_list.txt
    sort -k3,3nr ${RESULTS_DIR}/analysis_comp_microsat.txt.temp >> ${RESULTS_DIR}/microsatellite_cores_list.txt
    rm ${RESULTS_DIR}/analysis_comp_microsat.txt.temp

    printf "Comp_ID\tRepresentative\tMax abundance\n" >  ${RESULTS_DIR}/stretchAT_cores_list.txt
    sort -k3,3nr ${RESULTS_DIR}/analysis_comp_stretchAT.txt.temp >> ${RESULTS_DIR}/stretchAT_cores_list.txt
    rm ${RESULTS_DIR}/analysis_comp_stretchAT.txt.temp

    printf "Comp_ID\tRepresentative\tMax abundance\n" >  ${RESULTS_DIR}/other_cores_list.txt
    sort -k3,3nr ${RESULTS_DIR}/analysis_comp_others.txt.temp >> ${RESULTS_DIR}/other_cores_list.txt
    rm ${RESULTS_DIR}/analysis_comp_others.txt.temp

      end=`date +%s`
      elapsed=`expr $end - $begin`
      begin=`date +%s`
      echo "De novo classification time (in seconds): $elapsed \n"
fi

total_end=`date +%s`
total_elapsed=`expr $total_end - $start`
echo "Total time (in h:m:s): $(printf '%02d:%02d:%02d\n' $(($total_elapsed/3600)) $(($total_elapsed%3600/60)) $(($total_elapsed%60)))"
