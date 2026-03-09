#!/bin/bash
# TE analysis pipeline

begin=`date +%s`

#region LOAD THE ENVIRONMENT VARIABLES AND PATHS
if [[ -n "$1" && "$1" == *.sh ]]; then
    source "$1"
    ENV=$1
else
    source environment.sh
    ENV="environment.sh"
fi
#endregion

#region BUILD TE LIBRARY AND BINARIES
if [[ -z "${SKIP_TE_LIBRARIES}" ]]; then
  # Extraction des TEs de la base de données Dfam
  echo "Extraction of the TEs from Dfam ..."

  source ${VENV_FAMDB}
  #Parameters :
  #-i : path to the famdb installation
  #families : command to extract the families
  #--include-class-in-name : include the class of the TE in the name
  #--curated : only curated families
  #--descendants : include descendants of the specified species
  #--ancestors : include ancestors of the specified species
  #"${SPE_NAME}" : species name
  #--format fasta_name : output format with fasta header containing the TE name and description
  ${FAMDB_BIN} -i ${LIBRARY_DIR}/famdb/ families \
  --include-class-in-name \
  --curated \
  --descendants \
  --ancestors \
  "${SPE_NAME}" --format fasta_name \
   | sed 's/#/\t/g' \
   | sed 's/ @/\t/g' \
   | sed 's/^>/>dfam_/g' > ${DFAM_FA}

  #Build Bowtie genome (here that's the TEs)
  ${BOWTIE_BUILD} ${DFAM_FA} TE_Dfam_${SPE}


   sed 's/\t/_/g; s/ /_/g; s/\[/_/g; s/\]/_/g ; s/:/_/g; s/,/_/g ; s/\//_/g; s/(/_/g; s/)/_/g' \
   ${DFAM_FA} > ${DFAM_FA}.no_tab.fa

   ## TeTools to know the count of each TE
   awk '$1 ~ /^>/ {print $0, $2}' FS=['\t'_] ${DFAM_FA}.no_tab.fa | cut -c 2- > ${SPE_DIR}/rosette.fa

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Build TE library time (in seconds): $elapsed \n"
fi
#endregion

source ${WORK_DIR}/venv/bin/activate

#region BUILD THE BINARIES
if [[ -z "${SKIP_BUILD_RUST}" ]]; then
  ##Build the bin for the graph
  echo "Building the bin for the graph..."
  cargo build --release --manifest-path ${WORK_DIR}/graph/Cargo.toml
  cp ${WORK_DIR}/graph/target/release/graph ${WORK_DIR}/graph.exe

  ##Build the bin for the gene_finder_de_novo function
  echo "Building the bin for the gene_finder_de_novo function..."
  cargo build --release --manifest-path ${WORK_DIR}/gene_finder_de_novo/Cargo.toml
  cp ${WORK_DIR}/gene_finder_de_novo/target/release/gene_finder_de_novo ${WORK_DIR}/gene_finder_de_novo.exe

  end=`date +%s`
    elapsed=`expr $end - $begin`
    begin=`date +%s`
    echo "Build time (in seconds): $elapsed \n"
fi
#endregion

#region BOWTIE ALIGNMENT OF THE READS
if [[ -z "${SKIP_BOWTIE_READS}" ]]; then
  ##Align the read.fastp with bowtie2 to the Dfam TE library

  echo "Align the reads with bowtie2 to the Dfam TE library ..."


  mkdir -p ${RESULTS_DIR}/alignment/
  ####Bowtie2 alignment of the reads to the Dfam TE library
  ##Options used :
  # -a : report all alignments
  # -q : input files are in fastq format
  # -p 8 : use 8 threads
  # --very-sensitive : preset for very sensitive alignment
  # --dovetail : allow dovetailing alignments
  # -X 500 : maximum fragment length for valid paired-end alignments
  # -S : output SAM file
  # -1 : first read file
  # -2 : second read file paired with the first
  ${BOWTIE} --wrapper basic-0 \
      -a \
      -q \
      -p 8 \
      --time \
      --very-sensitive \
      -x TE_Dfam_${SPE} \
      --dovetail \
      -X 500 \
      -S ${RESULTS_DIR}/alignment/READS.sam \
      -1 ${READS_1}.fastp \
      -2 ${READS_2}.fastp &> ${RESULTS_DIR}/alignment/bowtie2_output_reads.txt

  #Convert SAM to BAM
  #Parameters :
  #-bS : input is SAM, output is BAM
  #-F 256 : filter out alignments with flag 256 (not primary alignments).
  # Bowtie2 randomly assigns primary alignments when multiple alignments have the same score
  # Keeping only primary alignments avoid counting multiple times the same read in downstream analyses
  ${SAMTOOLS_BIN} view -bS -F 256 \
      ${RESULTS_DIR}/alignment/READS.sam \
      > ${RESULTS_DIR}/alignment/READS.bam
  ${SAMTOOLS_BIN} sort \
      ${RESULTS_DIR}/alignment/READS.bam \
      -o ${RESULTS_DIR}/alignment/READS_sorted.bam
  ${SAMTOOLS_BIN} index \
      ${RESULTS_DIR}/alignment/READS_sorted.bam


  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Bowtie alignment of the reads time (in seconds): $elapsed \n"

fi
#endregion

#region BOWTIE ALIGNMENT OF THE COMPONENTS
if [[ -z "${SKIP_BOWTIE_COMPS}" ]]; then

  ##The number of components to compute
   MAXI=0 #$(ls ${BASE_DIR}/comp*.txt | wc -l)
      for file in ${BASE_DIR}/comp*.txt; do
        MAXI=$(( $MAXI + 1 ))
      done
      echo "Number of comps : ${MAXI}"


  # Print percentage every 1% (or every component if MAXI < 100)
    step=$(( MAXI / 100 ))
    if (( step == 0 )); then
        step=1
    fi

  ##Compute the analysis over every component
  for ((i=0; i<$MAXI; i++))
  do

      if (( i % step == 0 )); then
          PERCENTAGE=$(( (i * 100) / MAXI ))
          echo "Processing component $i of $MAXI ($PERCENTAGE% complete)"
      fi

      mkdir -p ${BASE_DIR}/alignment_${i}

      ##Convert the unitigs into reads for bowtie2

    python3 ${WORK_DIR}/reads_to_align.py \
        ${BASE_DIR}/comp${i}.txt \
        ${BASE_DIR}/comp${i}.fa \
        0
    ## -f : input files are in fasta format
      ${BOWTIE} --wrapper basic-0 \
          -a \
          -f \
          -p 1 \
          --time \
          --very-sensitive \
          -x TE_Dfam_${SPE} \
          --dovetail \
          -X 500 \
          -S ${BASE_DIR}/alignment_${i}/comp${i}.sam \
          -U ${BASE_DIR}/comp${i}.fa &> ${BASE_DIR}/alignment_${i}/bowtie2_output.txt


      ${SAMTOOLS_BIN} view -bS  \
          ${BASE_DIR}/alignment_${i}/comp${i}.sam \
          > ${BASE_DIR}/alignment_${i}/comp${i}.bam

      python3 ${WORK_DIR}/filter_bam.py \
          ${BASE_DIR}/alignment_${i}/comp${i}.bam \
          ${BASE_DIR}/aligned_Dfam${i}.txt

      :> ${BASE_DIR}/aligned_ref_AS${i}.txt
      :> ${BASE_DIR}/aligned_inter${i}.txt
      :> ${BASE_DIR}/aligned_Rb${i}.txt

      ##Annotate the components nodes with the TE and ref genes informations
      python3 ${WORK_DIR}/add_ref_TE.py \
              ${BASE_DIR}/comp${i}.txt \
              ${BASE_DIR}/aligned_ref_AS${i}.txt \
              ${BASE_DIR}/aligned_Dfam${i}.txt \
              ${BASE_DIR}/aligned_Rb${i}.txt \
              ${BASE_DIR}/aligned_inter${i}.txt \
              ${RESULTS_DIR}/graph/hc_1_hc_2_k${K}.abundance

  done

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Bowtie alignment of the reads time (in seconds): $elapsed \n"
fi
#endregion

#region ALIGNMENT OF THE UNITIGS
if [[ -z "${SKIP_BOWTIE_UNITIGS}" ]]; then
  mkdir -p ${BASE_DIR}/STAR_alignment_cons_dfam_all

  #Convert the unitigs into reads for star
  cp ${RESULTS_DIR}/graph/outputNodes.txt ${BASE_DIR}/all_unitigs.txt
  python3 ${WORK_DIR}/reads_to_align.py ${BASE_DIR}/all_unitigs.txt ${BASE_DIR}/all_unitigs.fa 0

  ###Bowtie2 alignment of the unitigs to the Dfam TE library
  ${BOWTIE} --wrapper basic-0 \
      -a \
      -f \
      -p 8 \
      --time \
      --very-sensitive \
      -x TE_Dfam_${SPE} \
      --dovetail \
      -X 500 -S ${BASE_DIR}/alignment_unitigs.sam \
      -U ${BASE_DIR}/all_unitigs.fa &> ${BASE_DIR}/bowtie2_output_unitigs.txt

  ${SAMTOOLS_BIN} view -bS \
      ${BASE_DIR}/alignment_unitigs.sam \
      > ${BASE_DIR}/alignment_unitigs.bam

    echo "Filter the reads to only keep the best mapped ones"

  python3 ${WORK_DIR}/filter_bam.py \
      ${BASE_DIR}/alignment_unitigs.bam \
      ${BASE_DIR}/aligned_Dfam_unitigs.txt

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Bowtie alignment of the unitigs time (in seconds): $elapsed \n"
fi
#endregion

#region ANALYSIS OF THE UNITIGS ALIGNMENT
if [[ -z "${SKIP_ANA_ALIGN}" ]]; then


#    awk '{print $1"\t"$20}' FS=["\t:"] \
#    ${BASE_DIR}/aligned_Dfam_unitigs.txt \
#    | sort -u -t '_' -k2 -n > ${BASE_DIR}/AS.txt

    :> ${BASE_DIR}/aligned_Rb_all.txt
    :> ${BASE_DIR}/aligned_ref_AS_all.txt
    :> ${BASE_DIR}/aligned_inter_all.txt

    echo "Annotate the nodes with the TE"
    ##Annotate the components nodes with the TE and ref genes informations
    python3 ${WORK_DIR}/add_ref_TE.py \
            ${BASE_DIR}/all_unitigs.txt \
            ${BASE_DIR}/aligned_ref_AS_all.txt \
            ${BASE_DIR}/aligned_Dfam_unitigs.txt \
            ${BASE_DIR}/aligned_Rb_all.txt \
            ${BASE_DIR}/aligned_inter_all.txt \
            ${RESULTS_DIR}/graph/hc_1_hc_2_k${K}.abundance

    end=`date +%s`
    elapsed=`expr $end - $begin`
    begin=`date +%s`
    echo "Analysis of the unitigs alignment time (in seconds): $elapsed \n"
fi
#endregion

#region ANALYSIS OF THE NEIGHBOURS OF THE COMPONENTS
if [[ -z "${SKIP_ANA_NEIGH}" ]]; then
    ##The number of components to compute
     MAXI=0 #$(ls ${BASE_DIR}/comp*.txt | wc -l)
        for file in ${BASE_DIR}/comp*.txt; do
          MAXI=$(( $MAXI + 1 ))
        done
    echo "Number of comps : ${MAXI}"


    mkdir -p ${BASE_DIR}/genes_of_comp
    ${WORK_DIR}/gene_finder_de_novo.exe \
          ${BASE_DIR}/comp \
          ${BASE_DIR}/all_unitigs_annotated.nodes \
          ${BASE_DIR}/../graph/hc_1_hc_2_k${K}_C0.05.edges \
          ${BASE_DIR}/../graph/hc_1_hc_2_k${K}.abundance \
          ${BASE_DIR}/genes_of_comp \
          ${K} \
          ${MAXI}

    end=`date +%s`
    elapsed=`expr $end - $begin`
    begin=`date +%s`
    echo "Analysis of the neighbours of the components time (in seconds): $elapsed \n"
fi
#endregion

#region TE COUNTING WITH TECOUNT OR FEATURECOUNTS+BEDTOOLS
if [[ -z "${SKIP_COUNT_EXPRESSED_TE}" ]]; then

#    OLD  COMMAND WITH TECOUNT :
#    echo "TE counting with TeCount ..."
#    Parameters :
#    -rosette : file with the exact TE consensus names and an asigned label of their family
#    -column : column number in the rosette file containing the TE family names (column 1 doesn"t work)
#    -TE_fasta : fasta file with the TE consensus sequences
#    -count : output file with the counts of each TE family
#    -bowtie2 : use bowtie2 alignment
#    -RNA : first reads file
#    -RNApair : second reads file
#    -insert : max insert size of the paired-end reads, here 500 bp to match the bowtie2 parameters in over steps
        ${TECOUNT} -rosette ${SPE_DIR}/rosette.fa  \
            -column 2 \
            -TE_fasta ${DFAM_FA}.no_tab.fa \
              -count ${SPE_DIR}/count_TE_TECOUNT.txt \
            -bowtie2 \
            -RNA ${READS_1} \
            -RNApair ${READS_2} \
            -insert 500

    echo "TE counting using bowtie2 all alignment of the reads to the Dfam TE library ..."

    #${RESULTS_DIR}/alignment/READS_sorted.bam : aligned reads to the Dfam TE library already sorted
    #Creat a SAF annotation file from the Dfam TE fasta file
    #SAF format :
    #GeneID  Chr Start   End Strand (with GeneID from 1 to end, Chr = TE name here, Start = 0, End = length of the TE, Strand = +)
    awk 'BEGIN {ORS=""} NR == 1 {print $1"\t"; next} /^>/ {print 0"\t"length(seq)"\t+\n"$1"\t"; seq=""} !/^>/ {seq=seq$0} END {print 0"\t"length(seq)"\t+"}' ${DFAM_FA} \
    | sed 's/>//g' \
    | awk 'BEGIN {print "GeneID\tChr\tStart\tEnd\tStrand"} {print NR"\t"$0}' \
    > ${DFAM_FA}.saf


    #Parameters :
      #-a : annotation file in SAF format
      #-o : output file with the counts of each TE family
      #-F SAF : format of the annotation file is SAF
      #--fracOverlap : minimum fraction of read length that must overlap a feature to be counted
      #-M : multi-mapping : REMOVED because we want only primary alignments
      #--fraction : fractional counting for multi-mapping reads : REMOVED (linked to -M)
      #-p : paired-end reads
      #-T : number of threads
    ${WORK_DIR}/featureCounts \
        -a ${DFAM_FA}.saf \
        -o ${SPE_DIR}/count_TE.txt \
        -F SAF \
        --fracOverlap 0.5 \
        -p \
        -T 8 \
        ${RESULTS_DIR}/alignment/READS_sorted.bam


    #Compute the TE coverage with bedtools
    echo "File conversion to BED"

    #Convert fasta to bed format
    #BED format :
    #chrom start end name (with chrom = TE name, start = 0, end = length of the TE, name = TE name)
    awk 'BEGIN {ORS=""} NR == 1 {print $1"\t"; next} /^>/ {print 0"\t"length(seq) "\n"$1"\t"; seq=""} !/^>/ {seq=seq$0} END {print 0"\t"length(seq)}' ${DFAM_FA} \
    | sed 's/>//g' \
    | awk '{print $0"\t"$1}' \
    | sort -k1,1 -k2,2n \
    > ${DFAM_FA}.bed

    #Convert the sorted bam to bed format
    ${BEDTOOLS_BIN}  bamtobed \
        -i ${RESULTS_DIR}/alignment/READS_sorted.bam | sort -k1,1 -k2,2n \
        > ${RESULTS_DIR}/alignment/READS.sorted.bed

    #Compute the coverage (-sorted option should work but here it raises an error)
    ${BEDTOOLS_BIN} coverage \
        -a ${DFAM_FA}.bed \
        -b ${RESULTS_DIR}/alignment/READS.sorted.bed \
        -nonamecheck \
        > ${RESULTS_DIR}/TE_coverage.txt

    #Reformat the files to kept only the numbers of counts and coverage per TE family
    awk 'NR>2 {print $2"\t"$7}' ${SPE_DIR}/count_TE.txt | sort -k1,1 > ${RESULTS_DIR}/TE_counts_final.txt
    awk '{print $1"\t"$3"\t"$8}' ${RESULTS_DIR}/TE_coverage.txt | sort -k1,1 > ${RESULTS_DIR}/TE_length_coverage_final.txt

    #Join the files; they have the same TE order
    join -t $'\t' ${RESULTS_DIR}/TE_counts_final.txt ${RESULTS_DIR}/TE_length_coverage_final.txt \
        | awk '{printf "%s\t%.3f\t%s\t%.2f\t%.2f\n", $1,$4,$2,76*$2/$3,76*$2/$3/98}'\
        | LC_ALL=C sort -t $'\t' -k4,4gr \
        | awk 'BEGIN {print "TE\tBreath_Coverage\tCount\tDepth_Coverage\tRPKM"} {print $0}' \
        > ${RESULTS_DIR}/TE_coverage_count_ab.txt

    #Keep only the line such counts > 0 and breadth coverage > 0.5
    awk 'NR == 0 || $2 > 0.5 && $3 > 0 && $4 > 1' ${RESULTS_DIR}/TE_coverage_count_ab.txt > ${RESULTS_DIR}/TE_coverage_count_ab_filtered.txt

    #Keep only the line such counts > 0
    awk 'NR == 0 || $3 > 0 ' ${RESULTS_DIR}/TE_coverage_count_ab.txt > ${RESULTS_DIR}/TE_coverage_count_ab_without_filtered.txt



    end=`date +%s`
    elapsed=`expr $end - $begin`
    begin=`date +%s`
    echo "TE counting time (in seconds): $elapsed \n"

fi
#endregion

#region TE VS COMPONENTS ANALYSIS
if [[ -z "${SKIP_TE_VS_COMPS}" ]]; then

    echo "TE vs Components analysis ..."

    #Getting TE in comps

    MAXI=0 #$(ls ${BASE_DIR}/comp*.txt | wc -l)
    for file in ${BASE_DIR}/comp*.txt; do
      MAXI=$(( $MAXI + 1 ))
    done
    echo "Number of comps : ${MAXI}"
    :> ${RESULTS_DIR}/temp.txt
    for ((i=0; i<$MAXI; i++))
    do
      awk '$6 != "*"{print $6}' FS="\t" ${BASE_DIR}/comp${i}_annotated.nodes \
      | sed 's/; /\n/g' \
      | sort -u  >> ${RESULTS_DIR}/temp.txt
    done
    sort -u ${RESULTS_DIR}/temp.txt > ${RESULTS_DIR}/TE_in_comps.txt
    rm ${RESULTS_DIR}/temp.txt

    #Looking for TE that are in the counts file but not in the comps
    awk 'NR==FNR{a[$1];next} !($1 in a){print $0}' FS="\t" \
    ${RESULTS_DIR}/TE_in_comps.txt \
    ${RESULTS_DIR}/TE_coverage_count_ab_filtered.txt \
    > ${RESULTS_DIR}/TE_not_in_comps.txt

    #And the opposite : TE in comps but not in counts file
    awk 'NR==FNR{a[$1];next} !($1 in a){print $0}' FS="\t" \
    ${RESULTS_DIR}/TE_coverage_count_ab_filtered.txt \
    ${RESULTS_DIR}/TE_in_comps.txt \
    > ${RESULTS_DIR}/TE_in_comps_not_in_counts.txt

    #Now check the stats for             ${BASE_DIR}/all_unitigs.txt
    awk '$6 != "*"{print $6}' FS="\t"  ${BASE_DIR}/all_unitigs_annotated.nodes \
    | sed 's/; /\n/g' \
    | sort -u  > ${RESULTS_DIR}/TE_in_all_unitigs.txt

    #TE in counts but not in all unitigs
    awk 'NR==FNR{a[$1];next} !($1 in a){print $0}' FS="\t" \
    ${RESULTS_DIR}/TE_in_all_unitigs.txt \
    ${RESULTS_DIR}/TE_coverage_count_ab_filtered.txt \
    > ${RESULTS_DIR}/TE_not_in_all_unitigs.txt

    #TE in all unitigs but not in counts
    awk 'NR==FNR{a[$1];next} !($1 in a){print $0}' FS="\t" \
    ${RESULTS_DIR}/TE_coverage_count_ab_filtered.txt \
    ${RESULTS_DIR}/TE_in_all_unitigs.txt \
    > ${RESULTS_DIR}/TE_in_all_unitigs_not_in_counts.txt

#TE in all unitigs but not in counts without restrictions
    awk 'NR==FNR{a[$1];next} !($1 in a){print $0}' FS="\t" \
    ${RESULTS_DIR}/TE_coverage_count_ab_without_filtered.txt \
    ${RESULTS_DIR}/TE_in_all_unitigs.txt \
    > ${RESULTS_DIR}/TE_in_all_unitigs_not_in_counts_without_filtering.txt



    end=`date +%s`
    elapsed=`expr $end - $begin`
    begin=`date +%s`
    echo "TE counting time (in seconds): $elapsed \n"
fi
#endregion

#region CONS_SEQ_ANALYSIS
if [[ -z "${SKIP_CONS_SEQ_ANA}" ]]; then

    echo "Analysis of comp consensus prediction..."

    #Getting the sequences from ${BASE_DIR}/analysis_comp_potential_TE.txt
    awk 'BEGIN{i=1} NR>1 {print ">SEQ_"i"_"$1"_"$3"\n"$2; i=i+1}' FS="\t" ${BASE_DIR}/analysis_comp_potential_TE.txt > ${BASE_DIR}/comp_potential_TE.fa
    sort -k1,1n ${BASE_DIR}/analysis_comp_potential_TE.txt > ${BASE_DIR}/analysis_comp_potential_TE.txt.sorted

    #Align the consensus sequences to the Dfam TE library with bowtie2
    ${BOWTIE} --wrapper basic-0 \
      -a \
      -f \
      -p 8 \
      --time \
      --very-sensitive \
      -x TE_Dfam_${SPE} \
      --dovetail \
      -X 500 -S ${BASE_DIR}/alignment_consensus.sam \
      -U ${BASE_DIR}/comp_potential_TE.fa &> ${BASE_DIR}/bowtie2_output_consensus.txt

    ${SAMTOOLS_BIN} view -bS \
      ${BASE_DIR}/alignment_consensus.sam \
      > ${BASE_DIR}/alignment_consensus.bam

#    ${SAMTOOLS_BIN} view \
#          ${BASE_DIR}/alignment_consensus.sam \
#          > ${BASE_DIR}/aligned_consensus.txt

    python3 ${WORK_DIR}/filter_bam.py \
      ${BASE_DIR}/alignment_consensus.bam \
      ${BASE_DIR}/aligned_consensus.txt
    awk 'NR>1 {print $0}' FS="\t" ${BASE_DIR}/analysis_comp_potential_TE.txt > ${BASE_DIR}/comp_potential_TE_unsorted.txt
    awk 'NR>1  {print $0}' FS="\t" ${BASE_DIR}/analysis_comp_microsat.txt  >> ${BASE_DIR}/comp_potential_TE_unsorted.txt
    awk 'NR>1  {print $0}' FS="\t" ${BASE_DIR}/analysis_comp_stretchA.txt >> ${BASE_DIR}/comp_potential_TE_unsorted.txt
    echo -e "Comp_ID \t Seq_consensus \t Max abundance" > ${BASE_DIR}/analysis_comp_all.txt
    sort -k3,3gr ${BASE_DIR}/comp_potential_TE_unsorted.txt | sort -u >> ${BASE_DIR}/analysis_comp_all.txt

    echo -e "Comp_ID \t Seq_consensus \t Max abundance" > ${BASE_DIR}/analysis_comp_all.txt.sorted
    sort -k1,1n ${BASE_DIR}/comp_potential_TE_unsorted.txt  | sort -u >> ${BASE_DIR}/analysis_comp_all.txt.sorted


    #Get the ROC curves for the consensus sequences
    python3 ${WORK_DIR}/seq_cons_roc_curve.py \
        ${BASE_DIR}/analysis_comp_potential_TE.txt \
        ${RESULTS_DIR}/TE_coverage_count_ab_filtered.txt \
        ${BASE_DIR}/aligned_consensus.txt \
        ${BASE_DIR}/comp \
        ${RESULTS_DIR}/output_roc_curve \
        ${BASE_DIR}/all_unitigs_annotated.nodes


    #Now idem with the 100 first comps
    echo "Analysis of first 100 comp consensus prediction..."
    awk 'NR>1 && $1 < 100 {print $0}' FS="\t" ${BASE_DIR}/analysis_comp_potential_TE.txt > ${BASE_DIR}/comp_potential_TE_100_unsorted.txt
    awk 'NR>1 && $1 < 100 {print $0}' FS="\t" ${BASE_DIR}/analysis_comp_microsat.txt  >> ${BASE_DIR}/comp_potential_TE_100_unsorted.txt
    awk 'NR>1 && $1 < 100 {print $0}' FS="\t" ${BASE_DIR}/analysis_comp_stretchA.txt >> ${BASE_DIR}/comp_potential_TE_100_unsorted.txt

    echo -e "Comp_ID \t Seq_consensus \t Max abundance" > ${BASE_DIR}/analysis_comp_potential_TE_100.txt
    sort -k3,3gr ${BASE_DIR}/comp_potential_TE_100_unsorted.txt >> ${BASE_DIR}/analysis_comp_potential_TE_100.txt

    echo -e "Comp_ID \t Seq_consensus \t Max abundance" > ${BASE_DIR}/analysis_comp_potential_TE_100.txt.sorted
    sort -k1,1n ${BASE_DIR}/comp_potential_TE_100_unsorted.txt >> ${BASE_DIR}/analysis_comp_potential_TE_100.txt.sorted

    rm ${BASE_DIR}/comp_potential_TE_100_unsorted.txt
    #Getting the sequences from ${BASE_DIR}/analysis_comp_potential_TE_100.txt
    awk 'BEGIN{i=1} NR>1 {print ">SEQ_"i"_"$1"_"$3"\n"$2; i=i+1}' FS="\t" ${BASE_DIR}/analysis_comp_potential_TE_100.txt > ${BASE_DIR}/comp_potential_TE_100.fa
    #Align the consensus sequences to the Dfam TE library with bowtie2
    ${BOWTIE} --wrapper basic-0 \
      -a \
      -f \
      -p 8 \
      --time \
      --very-sensitive \
      -x TE_Dfam_${SPE} \
      --dovetail \
      -X 500 -S ${BASE_DIR}/alignment_consensus_100.sam \
      -U ${BASE_DIR}/comp_potential_TE_100.fa &> ${BASE_DIR}/bowtie2_output_consensus_100.txt

    ${SAMTOOLS_BIN} view -bS \
      ${BASE_DIR}/alignment_consensus_100.sam \
      > ${BASE_DIR}/alignment_consensus_100.bam

    python3 ${WORK_DIR}/filter_bam.py \
      ${BASE_DIR}/alignment_consensus_100.bam \
      ${BASE_DIR}/aligned_consensus_100.txt

    #Get the ROC curves for the consensus sequences
    python3 ${WORK_DIR}/seq_cons_roc_curve.py \
        ${BASE_DIR}/analysis_comp_potential_TE_100.txt \
        ${RESULTS_DIR}/TE_coverage_count_ab_filtered.txt \
        ${BASE_DIR}/aligned_consensus_100.txt \
        ${BASE_DIR}/comp \
        ${RESULTS_DIR}/output_roc_curve_100

  #Look at the max extended degree of the missed TEs to understand if they are in small components or isolated nodes
  #This solution is quick and dirty, it should be optimised with a python script later
  awk 'NR>1 {print $1}' ${RESULTS_DIR}/TE_not_in_comps.txt > ${RESULTS_DIR}/missing_TE_names.txt
  : > ${RESULTS_DIR}/max_extended_deg_missed_TE.txt
  while read p;
   do
     grep "$p" ${BASE_DIR}/all_unitigs_annotated.nodes \
     | awk 'BEGIN {max=0} $4>max{max=$4} END {print max}' >> ${RESULTS_DIR}/max_extended_deg_missed_TE.txt
  done < ${RESULTS_DIR}/missing_TE_names.txt

  #paste Te_not_in_comps.txt and add the max extended degree
  sed -i '1s/^/Max_ext_deg\n/' ${RESULTS_DIR}/max_extended_deg_missed_TE.txt
  paste ${RESULTS_DIR}/TE_not_in_comps.txt ${RESULTS_DIR}/max_extended_deg_missed_TE.txt \
    > ${RESULTS_DIR}/TE_not_in_comps_with_max_extended_deg.txt



    end=`date +%s`
    elapsed=`expr $end - $begin`
    begin=`date +%s`
    echo "Consensus sequences analysis time (in seconds): $elapsed \n"

    python3 ${WORK_DIR}/unitigs_roc_curve.py \
            ${RESULTS_DIR}/TE_coverage_count_ab_filtered.txt \
            ${BASE_DIR}/all_unitigs_annotated.nodes \
            ${RESULTS_DIR}/output_roc_curve_unitigs
fi
#endregion