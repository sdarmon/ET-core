#The goal of this function is to analyse the core.txt file containing sequences in order to compute the
#number of poly(A) (represented as five consecutive A's) in the sequences; the ratio of microsatellites
#(repeated sequences of 1 to 6 nucleotides) in the sequences; and the annotated genes that are intersecting
# with the sequences.
import sys
import re
from collections import defaultdict
from numpy.ma.core import sum

Arg = sys.argv[:]

if len(Arg) not in [6]:
    print("Use : " + Arg[0] + " comp_prefix abundance_graph seq_consensus output_prefix connecting_unitigs")
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

def count_poly_CG(seq):
    compt = 0
    for i in range(len(seq) - 4):
        if seq[i:i+5] == 'CCCCC' or seq[i:i+5] == 'GGGGG':
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

#For each core, we will compute its degree in the induced subgraph
#and its maximal neighbor.
adjacent = [[] for _ in range(nb_comps)]
ids_save={}
total_deg = 0
with open(Arg[5],'r') as f :
    for line in f:
        if len(line) < 2:
            break
        L=line[:-1].split('\t')
        a=int(L[0])
        b=int(L[1])
        seq=L[2]
        adjacent[a].append(b)
        adjacent[b].append(a)
        total_deg+=2

mean_deg=total_deg/nb_comps

max_neighboor=[[] for _ in range(nb_comps)]
max_connection=[0 for _ in range(nb_comps)]
for i in range(nb_comps):
    nb_seen = {}
    for neigh in adjacent[i]:
        if neigh not in nb_seen:
            nb_seen[neigh] = 1
        else:
            nb_seen[neigh]+=1
    if nb_seen:
        max_nb = max(nb_seen.values())
        for neigh in nb_seen:
            if nb_seen[neigh] == max_nb:
                max_neighboor[i].append(neigh)
                max_connection[i]=max_nb

#Total output should look like :
#ID \t Representative \t Repeat type \t TE likelihood \t Max extended degree \t Max abundance \t Number of adjacent extended-t-cores \t Strongly connected to
#Where ID is the ID of the core
#Repeat type is A/T stretch; C/G stretch; Microsat (regex) XX%; others
#TE likelihood is -; +; ++; +++; depending on the repeat type other metrics
#Strongly connected to is the "ID (number of paths between them)"

#Keep only what is before the last / of Arg[4] for base name
l = Arg[4].split('/')
base_name = l[0]
for i in range(1, len(l)-1):
    base_name=base_name+"/"+l[i]


with open(Arg[4] + "_microsat.txt.temp", 'w') as f1:
    with open(base_name + "/extended_t_cores_summary.tsv", 'w') as f5:
        f5.write("Id \t Representative \t Repeat_type \t TE_Score \t Max_extended_degree \t Max_abundance \t Core_connectivity \t Best_neighbour\n")
        with open(Arg[4] + "_stretchAT.txt.temp", 'w') as f2:
            with open(Arg[4] + "_stretchCG.txt.temp", 'w') as f4:
                with open(Arg[4] + "_others.txt.temp", 'w') as f3:
                    for i in range(nb_comps):
                        #Reading the sequences from the core.txt file and compute the average number of poly(A) and the ratio of microsatellites
                        seqs = []
                        joined_seqs = ''
                        total_poly = 0
                        total_CG = 0
                        total_length = 0
                        ab_max = 0
                        score = 0
                        deg_max = 0
                        type_rep = ""
                        with open(Arg[1] + str(i) + ".txt", 'r') as f:
                            for line in f:
                                    if len(line) < 2:
                                        break
                                    L=line[:-1].split('\t')
                                    seqs.append(L[1])
                                    joined_seqs+= ' ' + L[1]
                                    total_poly += count_poly(L[1])
                                    total_CG += count_poly_CG(L[1])
                                    total_length += len(L[1])
                                    deg_max=max(deg_max, int(L[2]))
                                    ab_max= max(ab_max, abundance[int(L[0])])
                        #We count the overall ratio of microsat in the core
                        m, seq_m, r = count_microsat(joined_seqs)

                        #We count the ratio of unitigs of the core having at least 0.5 microsat
                        nb_microsat = 0
                        for seq in seqs:
                            m_seq, seq_m_seq, r_seq = count_microsat(seq)
                            if r_seq >= 0.5:
                                nb_microsat+=1
                        r_seqs= nb_microsat/len(seqs)

                        if total_poly / len(seqs) >= 0.75:
                            f2.write(f"{i}\t{seq_consensium[i]}\t{ab_max}\n")
                            type_rep="A/T stretch"
                            score+=1
                        elif total_CG / len(seqs) >= 0.75:
                            f4.write(f"{i}\t{seq_consensium[i]}\t{ab_max}\n")
                            type_rep="C/G stretch"
                            #score-=1
                        if r >= 0.5 or r_seqs >= 0.5:
                            f1.write(f"{i}\t{seq_consensium[i]}\t{ab_max}\t{seq_m}\t{r}\n")
                            texte="Microsat "+str(seq_m)+" ("+str(round(r*100))+"%)"
                            if type_rep=="":
                                type_rep= texte
                            else:
                                type_rep=type_rep+" and "+texte
                        if total_poly / len(seqs) < 0.75 and  total_CG / len(seqs) < 0.75 and (r < 0.5 and r_seqs < 0.5) :
                            f3.write(f"{i}\t{seq_consensium[i]}\t{ab_max}\n")
                            type_rep = "Potential TE"
                            score+=1
                        if len(adjacent[i]) > mean_deg :
                            score+=1

                        neigh_text="."
                        if max_connection[i] > mean_deg/2 :
                            score+=1
                            link= max_connection[i]
                            neigh = adjacent[i]
                            neigh_text=f"{max_neighboor[i][0]}:{100*max_connection[i]/len(adjacent[i]):.1f}%"

                        #Update the score
                        score_texte=str(score)


                        #Write in f5
                        #ID \t Representative \t Repeat type \t TE likelihood \t Max extended degree \t Max abundance \t Number of adjacent extended-t-cores \t Strongly connected to
                        #Where ID is the ID of the core
                        #Repeat type is A/T stretch; C/G stretch; Microsat (regex) XX%; others
                        #TE likelihood is -; +; ++; +++; depending on the repeat type other metrics
                        #Strongly connected to is the "ID (number of paths between them)"
                        f5.write(f"{i}\t{seq_consensium[i]}\t{type_rep}\t{score_texte}\t{deg_max}\t{ab_max}\t{len(adjacent[i])}\t{neigh_text}\n")