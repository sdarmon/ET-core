import numpy as np
import sys


Arg = sys.argv[:]

if len(Arg) not in [2]:
    print("Use : " + Arg[0] + "connected_edges.txt")
    exit()

#Reading the file
with open(Arg[1], 'r') as f:
    for line in f:
        L = line[:-1].split('\t')
        if len(L) < 2 :
            continue
        set_nodes = []
        for el in L:
            if el not in set_nodes:
                set_nodes.append(el)
        for i in range(len(set_nodes)):
            for j in range(i+1, len(set_nodes)):
                edge1 = set_nodes[i]
                edge2 = set_nodes[j]
                print(edge1+"\t"+edge2)
