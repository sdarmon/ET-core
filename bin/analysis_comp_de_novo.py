#The goal of this function is to analyse the comp.txt file containing sequences in order to compute the
#number of poly(A) (represented as five consecutive A's) in the sequences; the ratio of microsatellites
#(repeated sequences of 1 to 6 nucleotides) in the sequences; and the annotated genes that are intersecting
# with the sequences.
import sys
import bisect


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
    vu = [[0 for _ in range(encode("TTTTT")+1)] for _ in range(len(seq))]
    for i in range(6,(encode("TTTTT")+1)):
        s = decode(i)
        if s == '' or is_a_microsat(s):
            continue
        l = len(s)
        for j in range(len(seq) - l):
            if seq[j:j+l] == s and seq[j+l:j+2*l] == s:
                for k in range(j,j+2*l):
                    vu[k][i] = 1
                    i_rev = encode(rev_comp(s))
                    vu[k][i_rev] = 1
    #Finding the highest sum of vu[*][i] for i in range(6,(encode("TTTTT")+1))
    m = 0
    seq_m = ''
    for i in range(6,(encode("TTTTT")+1)):
        s = decode(i)
        if s == '' or is_a_microsat(s):
            continue
        nb_copy = sum([vu[j][i] for j in range(len(seq))]) / len(s)
        if nb_copy > m:
            m = nb_copy
            seq_m = s

    #Finding the positions covered by a microsatellite
    pos_vu = [0 for _ in range(len(seq))]
    for i in range(len(seq)):
        for j in range(6,(encode("TTTTT")+1)):
            if vu[i][j] == 1:
                pos_vu[i] = 1
    return (m, seq_m, sum(pos_vu)/len(seq))

#Define the sorted lists, ranked according the max ab value of each comp (by decreasing order) with bisect
microsat_comps = []
polyA_comps = []
potential_TE_comps = []

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

    if total_poly / len(seqs) >= 1:
        bisect.insort(polyA_comps, ( -ab_max, i))
    if r >= 0.2:
        bisect.insort(microsat_comps, ( -ab_max, i))
    if total_poly / len(seqs) < 1 and r < 0.2 :
        bisect.insort(potential_TE_comps, ( -ab_max, i))

#Writing the output files
with open(Arg[4] + "_microsat.txt", 'w') as f:
    f.write("Comp_ID\tSeq_consensus\tMax abundance\n")
    for el in microsat_comps:
        f.write(f"{el[1]}\t{seq_consensium[el[1]]}\t{-el[0]}\n")

with open(Arg[4] + "_stretchA.txt", 'w') as f:
    f.write("Comp_ID\tSeq_consensus\tMax abundance\n")
    for el in polyA_comps:
        f.write(f"{el[1]}\t{seq_consensium[el[1]]}\t{-el[0]}\n")

with open(Arg[4] + "_potential_TE.txt", 'w') as f:
    f.write("Comp_ID\tSeq_consensus\tMax abundance\n")
    for el in potential_TE_comps:
        f.write(f"{el[1]}\t{seq_consensium[el[1]]}\t{-el[0]}\n")
