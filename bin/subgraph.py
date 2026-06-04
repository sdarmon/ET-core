import sys

Arg = sys.argv[:]

if len(Arg) not in [6,7]:
    print("Use : " + Arg[0] + " core.txt graph.nodes graph.edges -o prefix_output [fast_cons]")
    exit()

if len(Arg) == 7:
    fast=True
    d_max = 1000
else:
    fast=False
    d_max = 20
core = []
max_ab = 0
id_max_ad = 0
comp_id = set()
with open(Arg[1],"r") as f:
    for line in f :
        L = line[:-1].split('\t')
        if len(L) < 2 :
            continue
        core.append(L)
        comp_id.add(int(L[0]))
        a = int(L[8])
        if a > max_ab:
            max_ab=a
            id_max_ad=int(L[0])

nodes = []
ab = []
seq = []
with open(Arg[2],"r") as f:
    for line in f :
        L=line[:-1].split("\t")
        if len(L) < 2 :
            continue
        ab.append(int(L[8]))
        seq.append(L[1])
        nodes.append((line[:-1],[]))

edges = []
edges_ways = {}
with open(Arg[3],"r") as f:
    for line in f :
        L = line[:-1].split('\t')
        if len(L) < 2 :
            continue

        nodes[int(L[0])][1].append(int(L[1]))
        edges_ways[(int(L[0]),int(L[1]))] = L[2]

id_a_voir = [(id_max_ad,d_max,True,True,'R','F')]

path = [id_max_ad]
left_path= []
right_path =[]
vu = set()

while id_a_voir != []:
    id, d, B1,B2, way1,way2 = id_a_voir[0]
    id_a_voir.pop(0)
    if id in vu:
        continue
    vu.add(id)
    if d == 0:
        continue
    ab_max_1 = 0
    ab_max_2 = 0
    id_max_1 = -1
    id_max_2 = -1
    n_way1 = ''
    n_way2= ''
    for id2 in nodes[id][1]:
        if ab[id2]>ab_max_1 and edges_ways[(id,id2)][0] == way1 and B1:
            ab_max_1=ab[id2]
            id_max_1=id2
            n_way1=edges_ways[(id,id2)][1]
    for id2 in nodes[id][1]:
        if ab[id2]>ab_max_2 and edges_ways[(id,id2)][0] == way2 and B2:
            ab_max_2=ab[id2]
            id_max_2=id2
            n_way2=edges_ways[(id,id2)][1]
    for id2 in nodes[id][1]:
        b1=False
        b2=False
        if id2 == id_max_1:
            b1 = True
            left_path.append(id2)
        if id2 == id_max_2:
            b2 = True
            right_path.append(id2)
        if id2 in comp_id: #Depth first in core
            id_a_voir= [(id2,d, b1,b2,n_way1,n_way2)] + id_a_voir
        elif (fast and (b1 or b2)):
            id_a_voir.append((id2,d-1,b1,b2,n_way1,n_way2))
        elif (not fast) and (b1 or b2): #Priorité aux chemin d'ab max
            id_a_voir= [(id2,d-1, b1,b2,n_way1,n_way2)] + id_a_voir
        elif not fast :
                id_a_voir.append((id2,d-1,b1,b2,n_way1,n_way2))

#Max ab sequence
def rc(s):
    S=''
    d = {}
    d['A']='T'
    d['C']='G'
    d['G']='C'
    d['T']='A'
    for n in s[::-1]:
        S=S+d[n]
    return S

seq_path = seq[id_max_ad]
seqs_in_the_path= [seq[id_max_ad]]
final_way='F'
for i in range(len(right_path)):
    if i == 0 :
        a = id_max_ad
    else:
        a = right_path[i-1]
    b = right_path[i]
    if b==a :
        continue
    if edges_ways[(a,b)]=="FF":
        seq_path=seq_path+seq[b][40:]
        final_way='F'
    elif edges_ways[(a,b)]=="RR":
        seq_path=seq[b][:-40]+seq_path
        final_way='R'
    elif edges_ways[(a,b)]=="FR":
        seq_path=seq[b][:-40]+rc(seq_path)
        final_way='R'
    else:
        seq_path=rc(seq_path)+seq[b][40:]
        final_way='F'

if final_way=='R':
    seq_path=rc(seq_path)

for i in range(len(left_path)):
    if i == 0 :
        a = id_max_ad
    else:
        a = left_path[i-1]
    b = left_path[i]
    if b==a :
        continue
    if edges_ways[(a,b)]=="FF":
        seq_path=seq_path+seq[b][40:]
    elif edges_ways[(a,b)]=="RR":
        seq_path=seq[b][:-40]+seq_path
    elif edges_ways[(a,b)]=="FR":
        seq_path=seq[b][:-40]+rc(seq_path)
    else:
        seq_path=rc(seq_path)+seq[b][40:]

print(seq_path)
if not fast:
    with open(Arg[5]+"_"+str(d_max)+"_.neigh","w") as f:
        for id in vu:
            if id in comp_id:
             f.write(nodes[id][0]+"\t1\n")
            else:
                f.write(nodes[id][0]+"\t0\n")


    with open(Arg[5]+"_"+str(d_max)+"_.neigh.edges","w") as f:
        for id in vu:
            for id2 in vu:
                if id < id2 and id2 in nodes[id][1] :
                    f.write(str(id)+"\t"+str(id2)+"\n")
