#The goal of this function is to analyse the comp.txt file containing sequences in order to compute the
#number of poly(A) (represented as five consecutive A's) in the sequences; the ratio of microsatellites
#(repeated sequences of 1 to 6 nucleotides) in the sequences; and the annotated genes that are intersecting
# with the sequences.
import sys
import gc
import re
from collections import defaultdict
from numpy.ma.core import sum

Arg = sys.argv[:]

if len(Arg) not in [5]:
    print("Use : " + Arg[0] + " comp_prefix abundance_graph seq_consensus output_prefix")
    exit()


#Read the abundance graph file to get the number of components
abundance=[]
with open(Arg[2], 'r') as f:
    for line in f:
        abundance.append(float(line[:-1]))

#Read the seq_consensium file
nb_comps=0
seq_consensium=[]
with open(Arg[3], 'r') as f:
    for line in f:
        nb_comps+=1
        seq_consensium.append(line[:-1])


#Function that count how many poly(A) a sequence has
def count_poly(seq):
    compt = 0
    for i in range(len(seq) - 4):
        if seq[i:i+5] == 'AAAAA' or seq[i:i+5] == 'TTTTT':
            compt += 1
    return compt


#Function that count the ratio of C/G in a sequence
def count_CG(seq):
    countCG = 0
    countAT = 0
    for el in seq:
        if el == 'C' or el == 'G':
            countCG += 1
        elif el == 'A' or el == 'T':
            countAT += 1
    return countCG / (countCG + countAT)

#Function that encodes a sequence in a binary format (blank=0, A=1, C=2, G=3, T=4)
def encode(seq):
    enc = {'A':1, 'C':2, 'G':3, 'T':4}
    s = 0
    for i in range(len(seq)):
        s = s*5 + enc[seq[i]]
    return s

#Function that decodes a sequence from a binary format to a nucleotide sequence
def decode(s):
    dec = {0:' ', 1:'A', 2:'C', 3:'G', 4:'T'}
    seq = ''
    while s > 0:
        if s%5 == 0:
            return ''
        seq = dec[s%5] + seq
        s = s//5
    return seq

#Function that checks if a sequence is a microsatellite
def is_a_microsat(seq):
    for i in range(1,len(seq)//2+1):
        if len(seq) % i != 0:
            continue
        if seq == seq[:i] * (len(seq) // i):
            return True
    return False

#Function that computes the reverse complement of a sequence
def rev_comp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join([comp[el] for el in seq[::-1]])

#Function that counts the number of microsatellites in a sequence
def count_microsat(seq):
    if not seq:
        return (0, "", 0.0)
    vu = [0 for _ in range(len(seq))]
    # Dictionary to keep track of total bases covered by each motif
    coverage = defaultdict(int)

    # Regex: Find 1 to 6 ACGT characters, repeated at least twice total.
    # '?' ensures it finds the shortest core motif (e.g., "A" instead of "AA")
    pattern = re.compile(r'([ACGT]{2,20}?)\1+')

    # finditer finds non-overlapping matches in a single fast C-level pass
    for match in pattern.finditer(seq):
        motif = match.group(1)      # The repeating unit (e.g., "AT")
        full_match = match.group(0) # The full tandem repeat (e.g., "ATATAT")

        #If the size of the full match is lesser than 8, we skip it
        if len(full_match) < 8:
            continue

        # Mark all positions covered by this tandem repeat as seen
        for i in range(match.start(), match.end()):
            vu[i] = 1

        # Add the length of the tandem repeat to this motif's total coverage
        coverage[motif] += len(full_match)

    if not coverage:
        return (0, "", 0.0)

    # Find the motif with the highest coverage
    seq_m = max(coverage, key=coverage.get)
    max_covered_bases = coverage[seq_m]

    # Calculate return values to match your original output format
    m = max_covered_bases // len(seq_m) # Total copies of this motif
    ratio = sum(vu) / len(seq) # Ratio of sequence covered by microsatellites

    return (m, seq_m, ratio)



with open(Arg[4] + "_microsat.txt.temp", 'w') as f1:
    with open(Arg[4] + "_stretchAT.txt.temp", 'w') as f2:
        with open(Arg[4] + "_others.txt.temp", 'w') as f3:
            for i in range(nb_comps):
                #Reading the sequences from the comp.txt file and compute the average number of poly(A) and the ratio of microsatellites
                seqs = []
                joined_seqs = ''
                total_poly = 0
                total_length = 0
                ab_max = 0
                with open(Arg[1] + str(i) + ".txt", 'r') as f:
                    for line in f:
                            if len(line) < 2:
                                break
                            L=line[:-1].split('\t')
                            seqs.append(L[1])
                            joined_seqs+= ' ' + L[1]
                            total_poly += count_poly(L[1])
                            total_length += len(L[1])
                            ab_max= max(ab_max, abundance[int(L[0])])
                    m, seq_m, r = count_microsat(joined_seqs)
                    if total_poly / len(seqs) >= 0.8:
                        f2.write(f"{i}\t{seq_consensium[i]}\t{ab_max}\n")
                    if r >= 0.2:
                        f1.write(f"{i}\t{seq_consensium[i]}\t{ab_max}\t{seq_m}\t{r}\n")
                    if total_poly / len(seqs) < 0.8 and r < 0.2 :
                        f3.write(f"{i}\t{seq_consensium[i]}\t{ab_max}\n")