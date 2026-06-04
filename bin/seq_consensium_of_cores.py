#The goal of this function is to analyse the comp.txt file containing sequences in
# order to compute the consensus sequence of the sequences. Since the component is
# connected, the consensus sequence is defined as the sequences of the path containing
# the heaviest node and 

import sys

Arg = sys.argv[:]

if len(Arg) not in [7]:
    print("Use : " + Arg[0] + " comp_prefix graph.edges graph.ab output k nb_comps")
    exit()

nb_comps= int(Arg[6])
k=int(Arg[5])
seq_cons_of_comps=[]

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return (''.join(complement[base] for base in seq[::-1]))

def cat(seq1, seq2, way2, k):
    #We suppose that seq1 has a forward orientation
    if way2 == 'F':
        s2 = seq2[:]
    else:
        s2 = reverse_complement(seq2)

    return (seq1[:] + s2[k-1:])


#Reading the graph.edges file to get the neighbours of each node
#edge format : node1 + "\t" + node2 + "\t" + type + "\n"
graph = {}
with open(Arg[2], 'r') as f:
    for line in f:
        L = line[:-1].split('\t')
        node1 = int(L[0])
        node2 = int(L[1])
        typeE = L[2][0:2]  # Remove the newline character at the end

        if node1 not in graph:
            graph[node1] = []
        graph[node1].append((node2,typeE))

ab_table=[]
with open(Arg[3], 'r') as f:
    for line in f:
        ab_table.append(float(line[:-1]))

for i in range(nb_comps):
    #Initialization of the dictionaries
    comp = {}
    max_node = -1
    max_ab = -1
    #Reading the comp_annotated.nodes file into a list
    #Line format : ind +"\t" + seq + "\t" + str(weight) + "\n"
    with open(Arg[1]+str(i)+".txt", 'r') as f:
        for line in f:
            L=line[:-1].split('\t')
            ind = int(L[0])
            seq = L[1]
            ab = ab_table[ind]
            comp[ind] = seq
            if ab > max_ab:
                max_ab = ab
                max_node = ind
    if max_node == -1:
        continue
    nodes_seen = set()
    nodes_seen.add(max_node)
    forward_path = []
    way = 'F'
    head = max_node
    # Start from the node with the maximum weight

    while True:
        neighbours = graph[head]
        next_node = -1
        best_ab = -1
        new_way = 'N'

        for (neighbour, edge_type) in neighbours:
            if neighbour in comp and neighbour not in nodes_seen and edge_type[0] == way:
                if ab_table[neighbour] > best_ab:
                    next_node = neighbour
                    best_ab = ab_table[neighbour]
                    new_way = edge_type[1]
        if next_node == -1:
            break
        forward_path.append((next_node, new_way))
        head = next_node
        nodes_seen.add(next_node)

    #Idem for the reverse path
    reverse_path = []
    way = 'R'
    head = max_node

    while True:
        neighbours = graph[head]
        next_node = -1
        best_ab = -1
        new_way = 'N'

        for (neighbour, edge_type) in neighbours:
            if neighbour in comp and neighbour not in nodes_seen and edge_type[0] == way:
                if ab_table[neighbour] > best_ab:
                    next_node = neighbour
                    best_ab = ab_table[neighbour]
                    new_way = edge_type[1]
        if next_node == -1:
            break

        reverse_path.append((next_node, new_way))
        head = next_node
        nodes_seen.add(next_node)

    #Now we have the forward and reverse paths, we can compute the consensus sequence

    #Compute the consensus sequence
    consensus_seq = comp[max_node]
    for node, way in forward_path:
        consensus_seq = cat(consensus_seq, comp[node],way, k)
    consensus_seq= reverse_complement(consensus_seq)
    for node, way in reverse_path:
        consensus_seq = cat(consensus_seq, comp[node],way, k)
    seq_cons_of_comps.append(consensus_seq)


#Output the consensus sequences of the components
with open(Arg[4], 'w') as f:
    for seq in seq_cons_of_comps:
        f.write(f"{seq}\n")
