import sys

Arg = sys.argv[:]

if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " comp.txt graph.nodes graph.edges")
    exit()

comp = []
comp_id = set()
with open(Arg[1],"r") as f:
    for line in f :
        L = line[:-1].split('\t')
        if len(L) < 2 :
            continue
        comp.append(L)
        comp_id.add(int(L[0]))

nodes = []
with open(Arg[2],"r") as f:
    for line in f :
        if len(L) < 2 :
            continue
        nodes.append((line[:-1],[]))

edges = []
with open(Arg[3],"r") as f:
    for line in f :
        L = line[:-1].split('\t')
        if len(L) < 2 :
            continue
        nodes[int(L[0])][1].append(int(L[1]))
        nodes[int(L[1])][1].append(int(L[0]))

d = 20
id_a_voir = [(int(L[0]),d) for L in comp]
vu = set()

while id_a_voir != []:
    id, d = id_a_voir[0]
    id_a_voir.pop(0)
    if id in vu:
        continue
    vu.add(id)
    if d == 0:
        continue
    for id2 in nodes[id][1]:
        id_a_voir.append((id2,d-1))

with open(Arg[1]+".neigh","w") as f:
    for id in vu:
        if id in comp_id:
         f.write(nodes[id][0]+"\t1\n")
        else:
            f.write(nodes[id][0]+"\t0\n")


with open(Arg[1]+".neigh.edges","w") as f:
    for id in vu:
        for id2 in vu:
            if id < id2 and id2 in nodes[id][1] :
                f.write(str(id)+"\t"+str(id2)+"\n")