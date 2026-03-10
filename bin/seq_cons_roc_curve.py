from matplotlib.colors import Normalize, LogNorm
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.collections as mcoll

def colorline(
        x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0),
        linewidth=3, alpha=1.0):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)

    ax = plt.gca()
    ax.add_collection(lc)

    return lc


def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

Arg = sys.argv[:]

if len(Arg) not in [7]:
    print("Use : " + Arg[0] + " analysis_comp_potential_TE.txt TE_coverage_count_ab_filtered.txt aligned_Dfam_consensus.txt comp_annotating_prefix output_roc_curve.png unitigs.nodes")
    exit()

#Reading TE_coverage_count_ab_filtered.txt to get the set of true TEs (positives)
true_TE_set = set()
with open(Arg[2], 'r') as f:
    i=0
    for line in f:
        if i==0:
            i=1
            continue
        L = line.split('\t')
        if len(L) < 2 :
            continue
        TE = L[0]
        true_TE_set.add(TE)

#Reading analysis_comp_potential_TE.txt to get the TEs sequences with their abundance
predicted_TE_list = []
with open(Arg[1], 'r') as f:
    i=0
    for line in f:
        if i==0:
            i=1
            continue
        L = line.split('\t')
        if len(L) < 2 :
            continue
        TE = L[1]
        abundance = float(L[2])
        extended_degree = int(L[0])
        predicted_TE_list.append((TE, abundance,extended_degree))


#predicted_TE_list is already sorted by abundance (from high to low)

#Read aligned_Dfam_consensus.txt to get the seq numbers that match to which TE
seq2TE = {}
with open(Arg[3], 'r') as f:
    for line in f:
        L = line.split('\t')
        if len(L) < 2 :
            continue
        index = int(L[0].split('_')[1])-1
        TE = L[2]
        if TE == "*":
            continue
        if index not in seq2TE:
            #creat an empty set
            seq2TE[index] = set()
        seq2TE[index].add(TE)

#Now compute a dynamic ROC curve points, keeping track of TP and FP in function of abundance threshold (from high to low)

TP = 0
TP_te = 0
FP = 0
FN = len(true_TE_set)
ab_threshold_list = []  #Start with a threshold higher than the highest abundance
TP_list = []
TP_te_list = []
FP_list = []
FN_list = []
TE_seen = set()
TE_seen_once = set()
multi_set=0
dic_comp2index = {}
for index in range(len(predicted_TE_list)):
    TE_seq, abundance, ex_deg = predicted_TE_list[index]
    #Get the TEs matching to this sequence
    TE_set = seq2TE[index] if index in seq2TE else set()
    dic_comp2index[ex_deg] = index

    if len(TE_set) > 1:
        multi_set+=1
    elif len(TE_set) ==1:
        TE_seen_once.add(list(TE_set)[0])
    if len(TE_set) > 0:
        TP+=1
    else :
        FP+=1

    for TE in TE_set:
        if TE in TE_seen:
            continue
        if TE in true_TE_set:
            FN -= 1
            TP_te+=1

        TE_seen.add(TE)

    #Update the lists
    ab_threshold_list.append(abundance)
    TP_list.append(TP)
    TP_te_list.append(TP_te)
    FP_list.append(FP)
    FN_list.append(FN)

print("Number of multi-TE matching sequences: "+str(multi_set)+"/"+str(len(predicted_TE_list)))
print("Number of unique TE matched: "+str(len(TE_seen_once))+"/"+str(len(TE_seen)))
#Now compute the dynamic precision-recall curve points
precision_list = []
recall_list = []
for i in range(len(TP_list)):
    TP = TP_list[i]
    TP_te = TP_te_list[i]
    FP = FP_list[i]
    FN = FN_list[i]
    precision = TP / (TP + FP) if (TP + FP) > 0 else 1.0
    recall = TP_te / (TP_te + FN) if (TP_te + FN) > 0 else 0.0
    precision_list.append(precision)
    recall_list.append(recall)

#Plot the ROC curve and color the dot corresponding to each abundance threshold
# ab_threshold_list doit exister et avoir la même longueur que recall_list/precision_list
c = np.array(ab_threshold_list)

# Choisir la normalisation : LogNorm() si large plage, sinon Normalize()
norm = LogNorm(vmin=c.min(), vmax=c.max()) if (c > 0).all() and (c.max() / c.min() > 100) else Normalize(vmin=c.min(), vmax=c.max())
fig, ax = plt.subplots(figsize=(6, 6))
# forcer les axes entre 0 et 1
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)
sc = ax.scatter(recall_list, precision_list, c=c, cmap='plasma', norm=norm, s=40, edgecolors='k')
ax.set_xlabel('Recall')
ax.set_ylabel('Precision')
ax.set_title('Precision-Recall Curve for TE Prediction by Ab Thresholding')
ax.grid(True)

# colorbar comme bloc-légende
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('abundance threshold')

plt.savefig(Arg[5]+"by_abundance.png", dpi=300)
plt.close()


#Now do it again with extended degree as threshold using this new file :
#Reading analysis_comp_potential_TE.txt.sorted to get the TEs sequences with their abundance
predicted_TE_list_ex = []
with open(Arg[1]+".sorted", 'r') as f:
    i=0
    for line in f:
        if i==0:
            i=1
            continue
        L = line.split('\t')
        if len(L) < 2 :
            continue
        TE = L[1]
        abundance = float(L[2])
        extended_degree = int(L[0])
        predicted_TE_list_ex.append((TE, abundance,extended_degree))

#Now compute a dynamic ROC curve points, keeping track of TP and FP in function of extended degree threshold (from high to low)
TP = 0
TP_te = 0
FP = 0
FN = len(true_TE_set)
ex_deg_threshold_list = []
TP_list = []
TP_te_list = []
FP_list = []
FN_list = []
TE_seen = set()
for index in range(len(predicted_TE_list_ex)):
    TE_seq, abundance, ex_deg = predicted_TE_list_ex[index]
    #Get the TEs matching to this sequence
    TE_set = seq2TE[dic_comp2index[ex_deg]] if dic_comp2index[ex_deg] in seq2TE else set()

    if len(TE_set) > 0:
        TP+=1
    else :
        FP+=1

    for TE in TE_set:
        if TE in TE_seen:
            continue
        if TE in true_TE_set:
            FN -= 1
            TP_te+=1

        TE_seen.add(TE)

    #Update the lists
    ex_deg_threshold_list.append(ex_deg)
    TP_list.append(TP)
    TP_te_list.append(TP_te)
    FP_list.append(FP)
    FN_list.append(FN)

#Now compute the dynamic precision-recall curve points
precision_list = []
recall_list = []
for i in range(len(TP_list)):
    TP = TP_list[i]
    TP_te = TP_te_list[i]
    FP = FP_list[i]
    FN = FN_list[i]
    precision = TP / (TP + FP) if (TP + FP) > 0 else 1.0
    recall = TP_te / (TP_te + FN) if (TP_te + FN) > 0 else 0.0
    precision_list.append(precision)
    recall_list.append(recall)

#Plot the ROC curve and color the dot corresponding to each extended degree threshold
# ex_deg_threshold_list doit exister et avoir la même longueur que recall_list/precision_list

c = np.array(ex_deg_threshold_list)
norm = LogNorm(vmin=c.min(), vmax=c.max()) if (c > 0).all() and (c.max() / c.min() > 100) else Normalize(vmin=c.min(), vmax=c.max())

fig, ax = plt.subplots(figsize=(6, 6))
# forcer les axes entre 0 et 1
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)
sc = ax.scatter(recall_list, precision_list, c=c, cmap='plasma', norm=norm, s=40, edgecolors='k')
ax.set_xlabel('Recall')
ax.set_ylabel('Precision')
ax.set_title('Precision-Recall Curve for TE Prediction by Extended Degree Thresholding')
ax.grid(True)

# forcer les axes entre 0 et 1
# colorbar attachée à la figure/sous-figure
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label('extended degree threshold')

fig.savefig(Arg[5] + "by_extended_degree.png", dpi=300)
plt.close(fig)


##Now we get the same curves but by taking all nodes from the comps
seq2TE={}
size=[]
size_dic = {}
ex_deg_2_actual_ex_deg = {}
nodes_seen = set()
min_deg= 10000
for i in range(len(predicted_TE_list)):
    TE_seq, abundance, ex_deg = predicted_TE_list[i]
    max_ed = 0
    with open(Arg[4]+str(ex_deg)+"_annotated.nodes", 'r') as f:
        size.append(0)
        for line in f:
            L = line.split('\t')
            if len(L) < 2 :
                continue
            size[-1]+=1
            nodes_seen.add(int(L[0]))
            index = i
            deg = int(L[3])
            if deg < min_deg:
                min_deg = deg
            if deg > max_ed:
                max_ed = deg
            TE = L[5]
            if TE == "*":
                continue
            if index not in seq2TE:
                #creat an empty set
                seq2TE[index] = set()
            for el in TE.split("; "):
                seq2TE[index].add(el)
    size_dic[ex_deg]=size[-1]
    ex_deg_2_actual_ex_deg[ex_deg] = max_ed


#Read the unitigs file
dic_additional_nodes = {}
dic_index2TE = {}
ex_deg_unitig_count = [0 for _ in range(min_deg)]
with (open(Arg[6], 'r') as f):
    for line in f:
        L = line.split('\t')
        if len(L) < 2 :
            continue
        index = int(L[0])
        ex_deg = int(L[3])
        TE = L[5]
        if index not in nodes_seen and ex_deg < min_deg:
            if ex_deg not in dic_additional_nodes:
                dic_additional_nodes[ex_deg] = []
            dic_additional_nodes[ex_deg].append(index)
            dic_index2TE[index] = TE
            ex_deg_unitig_count[ex_deg] += 1
#Now compute a dynamic ROC curve points, keeping track of TP and FP in function of extended degree threshold (from high to low)

TP = 0
TP_te = 0
FP = 0
FN = len(true_TE_set)
ab_threshold_list = []  #Start with a threshold higher than the highest abundance
TP_list = []
TP_te_list = []
FP_list = []
FN_list = []
TE_seen = set()

dic_comp2index = {}
for index in range(len(predicted_TE_list)):
    TE_seq, abundance, ex_deg = predicted_TE_list[index]
    #Get the TEs matching to this sequence
    TE_set = seq2TE[index] if index in seq2TE else set()
    dic_comp2index[ex_deg] = index

    if len(TE_set) > 0:
        TP+=1
    else :
        FP+=1

    for TE in TE_set:
        if TE in TE_seen:
            continue
        if TE in true_TE_set:
            FN -= 1
            TP_te+=1

        TE_seen.add(TE)

    #Update the lists
    ab_threshold_list.append(abundance)
    TP_list.append(TP)
    TP_te_list.append(TP_te)
    FP_list.append(FP)
    FN_list.append(FN)

#Now compute the dynamic precision-recall curve points
precision_list = []
recall_list = []
for i in range(len(TP_list)):
    TP = TP_list[i]
    TP_te = TP_te_list[i]
    FP = FP_list[i]
    FN = FN_list[i]
    precision = TP / (TP + FP) if (TP + FP) > 0 else 1.0
    recall = TP_te / (TP_te + FN) if (TP_te + FN) > 0 else 0.0
    precision_list.append(precision)
    recall_list.append(recall)

#Plot the ROC curve and color the dot corresponding to each abundance threshold
# ab_threshold_list doit exister et avoir la même longueur que recall_list/precision_list
c = np.array(ab_threshold_list)

# Choisir la normalisation : LogNorm() si large plage, sinon Normalize()
norm = LogNorm(vmin=c.min(), vmax=c.max()) if (c > 0).all() and (c.max() / c.min() > 100) else Normalize(vmin=c.min(), vmax=c.max())

#Normalise the size list into a np array from the size 40 to 240
size = np.array(size)
size = (size - size.min()) / (size.max() - size.min()) * 200 + 40


fig, ax = plt.subplots(figsize=(6, 6))
# forcer les axes entre 0 et 1
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)
sc = ax.scatter(recall_list, precision_list, c=c, cmap='plasma', norm=norm, s=size, edgecolors='k')
ax.set_xlabel('Recall')
ax.set_ylabel('Precision')


ax.set_title('Precision-Recall Curve for TE Prediction by Ab Thresholding and comps')
ax.grid(True)

# colorbar comme bloc-légende
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('abundance threshold')

plt.savefig(Arg[5]+"all_nodes_of_comp_by_abundance.png", dpi=300)
plt.close()

#Now do it again with extended degree as threshold using this new file :
TP = 0
TP_te = 0
FP = 0
FN = len(true_TE_set)
ex_deg_threshold_list = []
TP_list = []
TP_te_list = []
FP_list = []
FN_list = []
size_ex_deg=[]
TE_seen = set()
for index in range(len(predicted_TE_list_ex)):
    TE_seq, abundance, ex_deg = predicted_TE_list_ex[index]
    #Get the TEs matching to this sequence
    TE_set = seq2TE[dic_comp2index[ex_deg]] if dic_comp2index[ex_deg] in seq2TE else set()
    size_ex_deg.append(size_dic[ex_deg])
    if len(TE_set) > 0:
        TP+=1
    else :
        FP+=1

    for TE in TE_set:
        if TE in TE_seen:
            continue
        if TE in true_TE_set:
            FN -= 1
            TP_te+=1

        TE_seen.add(TE)

    #Update the lists
    ex_deg_threshold_list.append(ex_deg)
    TP_list.append(TP)
    TP_te_list.append(TP_te)
    FP_list.append(FP)
    FN_list.append(FN)

#convert the dic_additional_nodes keys to a sorted list from high to low
dic_keys = dic_additional_nodes.keys()
sorted_keys = sorted(dic_keys, reverse=True)
for deg in sorted_keys:
    for index in dic_additional_nodes[deg]:
        TE_set = dic_index2TE[index].split("; ") if dic_index2TE[index] != "*" else set()
        if len(TE_set) > 0:
            TP+=1
        else :
            FP+=1

        for TE in TE_set:
            if TE in TE_seen:
                continue
            if TE in true_TE_set:
                FN -= 1
                TP_te+=1

            TE_seen.add(TE)

        #Update the lists
        ab_threshold_list.append(abundance)
        TP_list.append(TP)
        TP_te_list.append(TP_te)
        FP_list.append(FP)
        FN_list.append(FN)


#Now compute the dynamic precision-recall curve points
precision_list = []
recall_list = []
for i in range(len(TP_list)):
    TP = TP_list[i]
    TP_te = TP_te_list[i]
    FP = FP_list[i]
    FN = FN_list[i]
    precision = TP / (TP + FP) if (TP + FP) > 0 else 1.0
    recall = TP_te / (TP_te + FN) if (TP_te + FN) > 0 else 0.0
    precision_list.append(precision)
    recall_list.append(recall)

#Plot the ROC curve and color the dot corresponding to each extended degree threshold
# ex_deg_threshold_list doit exister et avoir la même longueur que recall_list/precision_list
#complete ex_deg_threshold_list with the additional nodes extended degree (which is the key of dic_additional_nodes) and complete the size_ex_deg list
actual_ex_deg = [ex_deg_2_actual_ex_deg[ex_deg] for ex_deg in ex_deg_threshold_list]
for i in range(min_deg-1,-1,-1):
    actual_ex_deg+=[i]*ex_deg_unitig_count[i]
c = np.array(actual_ex_deg)
#Complete size_ex_deg with 1
size_ex_deg += [1]*(len(recall_list)-len(size_ex_deg))
#Normalise the size list into a np array from the size 40 to 240
size_ex_deg = np.array(size_ex_deg)
size_ex_deg = (size_ex_deg - size_ex_deg.min()) / (size_ex_deg.max() - size_ex_deg.min()) * 200 + 40

# Choisir la normalisation : LogNorm() si large plage, sinon Normalize()
norm = LogNorm(vmin=c.min(), vmax=c.max()) if (c > 0).all() and (c.max() / c.min() > 100) else Normalize(vmin=c.min(), vmax=c.max())
fig, ax = plt.subplots(figsize=(6, 6))
# forcer les axes entre 0 et 1
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)
sc = ax.scatter(recall_list, precision_list, c=c, cmap='plasma', norm=norm, s=size_ex_deg, edgecolors='k')
ax.set_xlabel('Recall')
ax.set_ylabel('Precision')


ax.set_title('Precision-Recall Curve for TE Prediction')
ax.grid(True)
# colorbar comme bloc-légende
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('Max extended degree')
plt.savefig(Arg[5]+"all_nodes_of_comp_by_extended_degree.png", dpi=300)
plt.close()
