##This function take a bam file and output a txt file where each sequence has only its
##maximal alignment scores in it. Thus, we are looking to keep only
##the best alignments (the ones that have an alignment score equal to the maximum for the read)
## for each sequence.

import pysam
import argparse
import os
import sys

def filter_bam(input_bam, output_bam):
    #open the input bam file
    input_bam_file = pysam.AlignmentFile(input_bam, "rb")
    #initialize a dictionary to store the tuples for each sequence:
    # (maximal alignment score, list of such line)
    max_scores = {}
    #iterate through the alignments in the input bam file
    for alignment in input_bam_file:
        #get the alignment score
        score = alignment.get_tag("AS")
        #get the read name
        read_name = alignment.query_name
        #if the read name is not in the dictionary, add it
        if read_name not in max_scores:
            max_scores[read_name] = (score, [alignment])
        #if the read name is in the dictionary, update the maximal score if needed
        else:
            if score > max_scores[read_name][0]:
                max_scores[read_name] = (score, [alignment])
            elif score == max_scores[read_name][0]:
                max_scores[read_name][1].append(alignment)
    #open the output bam file
    output_bam_file = pysam.AlignmentFile(output_bam, "wb", template=input_bam_file)
    #iterate through the dictionary and write the alignments to the output bam file
    for read_name in max_scores:
        for alignment in max_scores[read_name][1]:
            output_bam_file.write(alignment)
    #close the output bam file
    output_bam_file.close()
    #close the input bam file
    input_bam_file.close()

def filter_txt(input_bam, output_txt):
    #open the input bam file
    input_bam_file = pysam.AlignmentFile(input_bam, "rb")
    #initialize a dictionary to store the tuples for each sequence:
    # (maximal alignment score, list of such line)
    max_scores = {}
    #iterate through the alignments in the input bam file
    for alignment in input_bam_file:
        #get the alignment score
        score = alignment.mapping_quality
        #get the read name
        read_name = alignment.query_name
        #if the read name is not in the dictionary, add it
        if read_name not in max_scores:
            max_scores[read_name] = (score, [alignment])
        #if the read name is in the dictionary, update the maximal score if needed
        else:
            if score > max_scores[read_name][0]:
                max_scores[read_name] = (score, [alignment])
            elif score == max_scores[read_name][0]:
                max_scores[read_name][1].append(alignment)
    #close the input bam file
    input_bam_file.close()
    #open the output txt file where each line is the txt alignment
    output_txt_file = open(output_txt, "w")
    #iterate through the dictionary and write the alignments to the output bam file
    for read_name in max_scores:
        for alignment in max_scores[read_name][1]:
            output_txt_file.write(alignment.to_string() + "\n")
    #close the output bam file
    output_txt_file.close()

def main():
    #parse the arguments
    parser = argparse.ArgumentParser(description='Filter a bam file to keep only the best alignments for each sequence')
    parser.add_argument('input_bam', help='The input bam file')
    parser.add_argument('output_txt', help='The output txt file')
    #Add the option -b to put a flag as true for bam file
    parser.add_argument('-b', action='store_true', help='The input file is a bam file')
    args = parser.parse_args()
    #If an argument is missing print the help
    if args.input_bam is None or args.output_txt is None:
        parser.print_help()
        sys.exit(1)
    #check if the input bam file exists
    if not os.path.exists(args.input_bam):
        print("The input bam file does not exist")
        sys.exit(1)
    #call the filter_bam function if the flag is true , otherwise call the filter_txt function
    if args.b:
        filter_bam(args.input_bam, args.output_txt)
    else:
        filter_txt(args.input_bam, args.output_txt)

if __name__ == "__main__":
    main()
