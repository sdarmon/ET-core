#This code read a a comp%i.txt file and a seq_intersectionTE%i.txt file
#and edit the comp%i.txt file to add the TE reference in the seq_intersectionTE%i.txt file for every line
# that is present in the seq_intersectionTE%i.txt file


import sys

Arg = sys.argv[:]

if len(Arg) not in [7,8]:
    print("Use : " + Arg[0] + " comp%i.txt seq_intersectionRef%i.txt seq_intersectionConsensusDfam%i.txt seq_intersectionConsensusRb%i.txt seq_intergenes$i.txt file.abundance [threshold]")
    exit()

if len(Arg) == 8:
    threshold = int(Arg[7])
else:
    threshold = 0
#Reading intersectionRef%i.txt into a list to get the gene intersection with the comps
dic_index2gene = {}
exon_non_codant = {}
potential_intron = {}
"""
gene_biotype = ""
for info in L[14].split(';'):
    if info.startswith(" gene_biotype"):
        gene_biotype = info.split('"')[1]
        break
                
                
elif gene_biotype == "lncRNA" :
typed = "lncRNA"
"""


prop_intron = 0.35 #if a sequence has less than 30% of its length inside a CDS
# UTR, NCE, it is considered as an intron


with open(Arg[2], 'r') as f:
    for line in f:
        L = line.split('\t')
        len_seq = len(L[1])
        gene = L[20].split('"')[1] + L[5]#Old L[14] ; L[5] is "+" or "-"
        index = int(L[3].split('_')[1])
        position = L[14] #Old L[8]
        overlap = int(L[21]) #Old L[15]
        AS = int(L[22]) #Old L[16]
        if position in ["CDS","start_codon","stop_codon"] and overlap > prop_intron * len_seq :
            typed="CDS"
        elif position in ["five_prime_utr","three_prime_utr"] and overlap > prop_intron * len_seq:
            typed="UTR"
        elif position in ["exon"] and overlap > prop_intron * len_seq:
            #Case of a potential exon that is not in the CDS
            if index not in exon_non_codant:
                exon_non_codant[index] = [(gene,AS,overlap)]
            elif gene not in exon_non_codant[index]:
                exon_non_codant[index].append((gene,AS,overlap))
            continue
        elif position in ["gene", "transcript"]:
            if index not in potential_intron:
                potential_intron[index] = [(gene,AS,overlap)]
            elif gene not in potential_intron[index]:
                potential_intron[index].append((gene,AS,overlap))
            continue
        else: #Should not happen
            continue
        if index not in dic_index2gene:
            dic_index2gene[index] = [(gene,typed,AS,overlap)]
        elif (gene,typed,AS) not in dic_index2gene[index]:
            dic_index2gene[index].append((gene,typed,AS,overlap))
    #Now loop over the potential exon to add the exon to the gene list if it is not already in the list
    for index in exon_non_codant:
        for (gene,AS,overlap) in exon_non_codant[index]:
            if index not in dic_index2gene :
                dic_index2gene[index] = [(gene,"NCE",AS,overlap)]
            elif gene not in [el[0] for el in dic_index2gene[index]] :
                dic_index2gene[index].append((gene,"NCE",AS,overlap))


    #Now loop over the potential intron to add the intron to the gene list if it is not already in the list
    for index in potential_intron:
        for (gene,AS,overlap) in potential_intron[index]:
            if index not in dic_index2gene :
                dic_index2gene[index] = [(gene,"intron",AS,overlap)]
            elif gene not in [el[0] for el in dic_index2gene[index]] :
                dic_index2gene[index].append((gene,"intron",AS,overlap))



#Reading intersectionConsensus%i.txt into a list to get the TE consensus intersection with the comps
dic_index2consDfam = {}

with open(Arg[3], 'r') as f:
    for line in f:
        L = line.split('\t')
        gene = L[2]
        if gene == "*":
            continue
        index = int(L[0].split('_')[1]) #Here the index is the line of comp%i.txt file
        AS = int(L[13].split(':')[2])
        if index not in dic_index2consDfam:
            dic_index2consDfam[index] = [(gene,AS)]
        elif gene not in dic_index2consDfam[index]:
            dic_index2consDfam[index].append((gene,AS))


#Reading intersectionConsensusRB%i.txt into a list to get the TE consensus intersection with the comps
dic_index2consRb = {}

with open(Arg[4], 'r') as f:
    for line in f:
        L = line.split('\t')
        if len(L) < 2 :
            continue
        gene = L[2]
        if gene == "*":
            continue
        index = int(L[0].split('_')[1])
        if index not in dic_index2consRb:
            dic_index2consRb[index] = [gene]
        elif gene not in dic_index2consRb[index]:
            dic_index2consRb[index].append(gene)


#Reading intersectionGenes%i.txt into a list to get the gene intersection with the comps$i
dic_index2intergene = {}
with open(Arg[5], 'r') as f:
    for line in f:
        L = line.split('\t')
        if len(L) < 2 :
            continue
        index = int(L[0].split('_')[1])
        chromosome = L[2]
        deb = L[3]
        if chromosome == "*" :
            continue
        if index not in dic_index2intergene:
            dic_index2intergene[index] = [(chromosome,deb)]
        elif (chromosome,deb) not in dic_index2intergene[index]:
            dic_index2intergene[index].append((chromosome,deb))

#Reading the abundance file
abundance = []
with open(Arg[6], 'r') as f:
    for line in f:
        abundance.append(line.split('.')[0])


#Reading every line of comp%i.txt and adding the TE reference if the index of the line is in dic_index2TE
#Writing the new line in a new file comp%i.txt
#Here there is a point to be careful, the index of the line in comp%i.txt (a subset from the total number of nodes)
# is not the same as the index in seq_intersectionTE%i.txt (from 0 to the number of nodes of the comp)

with open(Arg[1], 'r') as f:
    with open(Arg[1].split('.')[0] + "_annotated.nodes", 'w') as f_out:
        i = 0 #Here is the index of the line in comp%i.txt
        for line in f:
            L = line[:-1].split('\t')
            ind=L[0]
            seq=L[1]
            node = int(L[0]) #Here is the actual index of the unitig of whole graph
            AS=0
            AS_dfam=0
            d=0
            str_gene = "*"
            str_consDfam = "*"
            str_consRb = "*"
            str_intergene = "*"

            if i in dic_index2gene:
                str_gene = "; ".join([el[0] + "@" + el[1] + "%" + str(el[3]) for el in dic_index2gene[i]])
                AS = max(el[2] for el in dic_index2gene[i])
                d+=4

            if i in dic_index2consDfam:
                d+=2
                str_consDfam = "; ".join([el[0] for el in dic_index2consDfam[i]])
                AS_dfam = max(el[1] for el in dic_index2consDfam[i])

            if i  in dic_index2consRb:
                str_consRb = "; ".join(dic_index2consRb[i])

            if i in dic_index2intergene:
                temp = [el[0] + ":" + el[1] for el in dic_index2intergene[i]]
                str_intergene = "; ".join(temp)

            if len(L) == 3:
                weight= L[2]
            else:
                weight = L[3]
            if int(weight) >= threshold:
                d+=1
            f_out.write(ind +"\t" + seq + "\t" + str(d) + "\t" + weight + "\t" + str_gene + "\t" + str_consDfam  + "\t" + str_consRb  + "\t" + str_intergene  + "\t" +  abundance[int(L[0])] + "\t" + str(AS) + "\t" + str(AS_dfam) + "\n")

            i += 1
