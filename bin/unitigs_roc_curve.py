from fontTools.misc.py23 import range
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

if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " TE_coverage_count_ab_filtered.txt unitig_annotated.nodes output_roc_curve.png")
    exit()

#Reading TE_coverage_count_ab_filtered.txt to get the set of true TEs (positives)
true_TE_set = set()
with open(Arg[1], 'r') as f:
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

#Reading unitig_annotated.nodes
dic_ext_deg_to_ids = {}
TEs=[]
sum_element_count = 0
with open(Arg[2], 'r') as f:
    for line in f:
        L = line.split('\t')
        if len(L) < 2 :
            continue
        id = int(L[0])
        extended_degree = int(L[3])
        if extended_degree >= 10:
            sum_element_count+=1
        if extended_degree not in dic_ext_deg_to_ids:
            dic_ext_deg_to_ids[extended_degree] = set()
        dic_ext_deg_to_ids[extended_degree].add(id)
        TEs.append(L[5])



#Now compute a dynamic ROC curve points, keeping track of TP and FP

TP_total = 0
FN = len(true_TE_set)
extended_degree_keys = sorted(dic_ext_deg_to_ids.keys(), reverse=True)  # Sort extended degrees from high to low
TP_list = []
TP_total_list = []
FP_list = []
FN_list = []
TE_seen = set()
TE_seen_once = set()
multi_set=0
dic_comp2index = {}
for deg in extended_degree_keys:
    id_set = dic_ext_deg_to_ids[deg]
    TE_set = set()
    TP = 0
    FP = 0
    for id in id_set:
        TE = TEs[id]
        if TE != "*":
            TP+=1
            for el in TE.split("; "):
                if el not in TE_seen:
                    TE_set.add(el)
        else:
            FP+=1

    if len(TE_set) > 1:
        multi_set+=1
    elif len(TE_set) ==1:
        TE_seen_once.add(list(TE_set)[0])

    for TE in TE_set:
        if TE in TE_seen:
            continue
        if TE in true_TE_set:
            TP_total+=1
            FN -= 1

        TE_seen.add(TE)

    #Update the lists for the ROC curve
    TP_list.append(TP)
    TP_total_list.append(TP_total)
    FP_list.append(FP)
    FN_list.append(FN)


#Now compute the dynamic precision-recall curve points
precision_list = []
recall_list = []
max_pre = 0
min_deg = 0
for i in range(len(TP_list)):
    TP = TP_list[i]
    TP_total = TP_total_list[i]
    FP = FP_list[i]
    FN = FN_list[i]
    precision = TP / (TP + FP) if (TP + FP) > 0 else 1.0
    recall = TP_total / (TP_total + FN) if (TP_total + FN) > 0 else 0.0
    if recall > 0.2 and precision > max_pre:
        max_pre = precision
        min_deg = extended_degree_keys[i]
    precision_list.append(precision)
    recall_list.append(recall)

#Plot the ROC curve and color the dot corresponding to each extended degree threshold by the value of the extended degree threshold
c = np.array(extended_degree_keys)

# Choisir la normalisation : LogNorm() si large plage, sinon Normalize()
norm = LogNorm(vmin=c.min(), vmax=c.max()) if (c > 0).all() and (c.max() / c.min() > 100) else Normalize(vmin=c.min(), vmax=c.max())
fig, ax = plt.subplots(figsize=(6, 6))
# forcer les axes entre 0 et 1
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)
sc = ax.scatter(recall_list, precision_list, c=c, cmap='plasma', norm=norm, s=40, edgecolors='k')

#Add the max_pre point in red with its extended degree threshold as label
ax.scatter(recall_list[precision_list.index(max_pre)], max_pre, c='green', s=100, edgecolors='k', label=f'Max Precision: {max_pre:.2f} at Extended Degree: {min_deg}')

ax.set_xlabel('Recall (by extended degree thresholding)')
ax.set_ylabel('Precision (by extended degree)')
ax.set_title('Precision-Recall Curve for TE Prediction')
ax.grid(True)

# colorbar comme bloc-légende
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('Extended degree threshold')


#Print the label of the max_pre point on top of the dot
ax.legend(loc='upper right')

fig.savefig(Arg[3] + "_by_extended_degree.png", dpi=300)
plt.close(fig)


##Now do the same plot but with interval of 100, divided by 10th-tiles
interval = 100
max_deg = max(extended_degree_keys)

interval_keys = [0]
corresponding_intervals = {}
corresponding_intervals[max_deg]=0
deg_table=[max_deg]
total_number_of_elements = sum_element_count
sum_elements = 0
for i in range(len(extended_degree_keys)):
    if extended_degree_keys[i] < 10 :
        interval_keys.append(len(dic_ext_deg_to_ids[extended_degree_keys[i]]))
        corresponding_intervals[extended_degree_keys[i]]=len(interval_keys)-1
        deg_table.append(extended_degree_keys[i])
        continue
    nb_elements = len(dic_ext_deg_to_ids[extended_degree_keys[i]])
    sum_elements += nb_elements
    while sum_elements > (total_number_of_elements/interval)*len(interval_keys) :
        interval_keys.append(0)
        deg_table.append(extended_degree_keys[i])
    interval_keys[-1] += nb_elements
    corresponding_intervals[extended_degree_keys[i]]=len(interval_keys)-1
    deg_table[-1] = extended_degree_keys[i]
    

TP_interval = [0] * len(interval_keys)
TP_total_interval = [0] * len(interval_keys)
FP_interval = [0] * len(interval_keys)
FN_interval = [len(true_TE_set)] * len(interval_keys)
for i in range(len(extended_degree_keys)):
    deg = extended_degree_keys[i]
    j = corresponding_intervals[deg]
    TP = TP_list[i]
    TP_total = TP_total_list[i]
    FP = FP_list[i]
    FN = FN_list[i]

    TP_interval[j] += TP
    FP_interval[j] += FP
    FN_interval[j] = FN
    TP_total_interval[j] = TP_total

precision_interval = []
recall_interval = []
for i in range(len(interval_keys)):
    TP = TP_interval[i]
    TP_total = TP_total_interval[i]
    FP = FP_interval[i]
    FN = FN_interval[i]
    precision = TP / (TP + FP) if (TP + FP) > 0 else 1.0
    recall = TP_total / (TP_total + FN) if (TP_total + FN) > 0 else 0.0
    precision_interval.append(precision)
    recall_interval.append(recall)


#Plot the ROC curve and color the dot corresponding to each extended degree threshold by the value of the extended degree threshold
c = np.array(deg_table)

# Choisir la normalisation : LogNorm() si large plage, sinon Normalize()
norm = LogNorm(vmin=c.min(), vmax=c.max()) if (c > 0).all() and (c.max() / c.min() > 100) else Normalize(vmin=c.min(), vmax=c.max())
fig, ax = plt.subplots(figsize=(6, 6))
# forcer les axes entre 0 et 1
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)
sc = ax.scatter(recall_interval, precision_interval, c=c, cmap='plasma', norm=norm, s=40, edgecolors='k')


ax.set_xlabel('Recall (by extended degree thresholding)')
ax.set_ylabel('Precision (by extended degree)')
ax.set_title('Precision-Recall Curve for TE Prediction')
ax.grid(True)

ax.scatter(recall_list[precision_list.index(max_pre)], max_pre, c='green', s=100, edgecolors='k', label=f'Max Precision: {max_pre:.2f} at Extended Degree: {min_deg}')

# colorbar comme bloc-légende
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('Extended degree threshold')

#Print the label of the max_pre point on top of the dot
ax.legend(loc='upper right')

fig.savefig(Arg[3] + "_by_extended_degree_interval_"+str(interval)+".png", dpi=300)
plt.close(fig)
