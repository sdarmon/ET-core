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


def size_to_proba(size, rate):
    # From a size number of pick elements (of probability rate) compute the probability to pick at least one of them
    return 1 - (1 - rate) ** size


def _parse_te_field(field):
    """Return a set of TE names from a core label; empty set = false positive."""
    if field is None:
        return set()
    if isinstance(field, (list, tuple, set)):
        return {str(x).strip() for x in field if str(x).strip() and str(x).strip() != "*"}
    s = str(field).strip()
    if s == "" or s == "*":
        return set()
    # allow "TE1; TE2" or "TE1;TE2"
    parts = [p.strip() for p in s.replace(";", "; ").split("; ")]
    return {p for p in parts if p and p != "*"}


def cumulative_pr_from_core_list(core_te_labels, true_te_set, unitig_te_labels=None):
    """
    Walk cores (then optional unitigs) in order and build cumulative precision / recall.

    Parameters
    ----------
    core_te_labels : list
        Ordered list, one entry per core. Examples:
          ["TE1", "", "TE2"]
          ["TE1; TE2", "", "TE3"]
    true_te_set : set of str
        Expressed TE consensus names (ground truth for recall).
    unitig_te_labels : list, optional
        Extra points after cores (e.g. unitigs with extended degree < t),
        same label format as cores.

    Returns
    -------
    precision_list, recall_list, TP_list, FP_list, TP_te_list, FN_list, n_cores
        n_cores is the number of core points (before unitig extension).
    """
    true_te_set = set(true_te_set)
    TP = 0
    TP_te = 0
    FP = 0
    FN = len(true_te_set)
    TE_seen = set()

    TP_list, FP_list, TP_te_list, FN_list = [], [], [], []
    precision_list, recall_list = [], []

    labels = list(core_te_labels)
    n_cores = len(labels)
    if unitig_te_labels is not None:
        labels = labels + list(unitig_te_labels)

    for field in labels:
        TE_set = _parse_te_field(field)
        if len(TE_set) > 0:
            TP += 1
        else:
            FP += 1

        for TE in TE_set:
            if TE in TE_seen:
                continue
            if TE in true_te_set:
                FN -= 1
                TP_te += 1
            TE_seen.add(TE)

        TP_list.append(TP)
        FP_list.append(FP)
        TP_te_list.append(TP_te)
        FN_list.append(FN)
        precision = TP / (TP + FP) if (TP + FP) > 0 else 1.0
        recall = TP_te / (TP_te + FN) if (TP_te + FN) > 0 else 0.0
        precision_list.append(precision)
        recall_list.append(recall)

    return precision_list, recall_list, TP_list, FP_list, TP_te_list, FN_list, n_cores


def plot_pr_from_core_list(
        core_te_labels,
        true_te_set,
        output_png="pr_from_core_list.png",
        ranks=None,
        title="Precision-Recall Curve for TE Prediction",
        rate=None,
        core_sizes=None,
        unitig_te_labels=None,
        unitig_ranks=None,
):
    """
    Plot a precision-recall curve from an ordered core?TE label list,
    optionally extended with unitigs (size 1 on the plot).

    ranks : optional list of scores to colour core points (e.g. max extended degree).
    unitig_te_labels / unitig_ranks : optional extension after cores (Ch.4-style U^{<t}).
    rate / core_sizes : optional null-model overlay (same spirit as seq_cons_roc_curve.py).
    """
    n_cores = len(core_te_labels)
    n_unitigs = len(unitig_te_labels) if unitig_te_labels is not None else 0

    if ranks is not None and len(ranks) != n_cores:
        raise ValueError(
            f"ranks length ({len(ranks)}) must match core_te_labels ({n_cores})"
        )
    if core_sizes is not None and len(core_sizes) != n_cores:
        raise ValueError(
            f"core_sizes length ({len(core_sizes)}) must match core_te_labels ({n_cores})"
        )
    if unitig_te_labels is not None:
        if unitig_ranks is None:
            raise ValueError("unitig_ranks is required when unitig_te_labels is set")
        if len(unitig_ranks) != n_unitigs:
            raise ValueError(
                f"unitig_ranks length ({len(unitig_ranks)}) must match "
                f"unitig_te_labels ({n_unitigs})"
            )

    precision_list, recall_list, TP_list, FP_list, TP_te_list, FN_list, n_cores = (
        cumulative_pr_from_core_list(
            core_te_labels, true_te_set, unitig_te_labels=unitig_te_labels
        )
    )

    n = n_cores + n_unitigs
    if ranks is None:
        ranks = list(range(n_cores, 0, -1))
    ranks_all = list(ranks)
    if unitig_ranks is not None:
        ranks_all = ranks_all + list(unitig_ranks)
    c = np.array(ranks_all, dtype=float)

    if core_sizes is None:
        size_cores = np.full(n_cores, 80.0)
    else:
        size_cores = np.array(core_sizes, dtype=float)
        if n_cores > 0 and size_cores.max() > size_cores.min():
            size_cores = (size_cores - size_cores.min()) / (size_cores.max() - size_cores.min()) * 200 + 40
        else:
            size_cores = np.full(n_cores, 80.0)
    # Unitigs: fixed marker size 1 (then scaled like a tiny point for visibility)
    size_unitigs = np.full(n_unitigs, 1.0)
    if n_unitigs > 0:
        # map size value 1 onto the same visual scale as cores (min bubble ~40)
        if core_sizes is not None and n_cores > 0:
            raw = np.concatenate([np.array(core_sizes, dtype=float), size_unitigs])
            if raw.max() > raw.min():
                size = (raw - raw.min()) / (raw.max() - raw.min()) * 200 + 40
            else:
                size = np.concatenate([size_cores, np.full(n_unitigs, 40.0)])
        else:
            size = np.concatenate([size_cores, np.full(n_unitigs, 40.0)])
    else:
        size = size_cores

    last_r_cores = recall_list[n_cores - 1] if n_cores > 0 else 0.0
    last_r = recall_list[-1] if recall_list else 0.0
    print(
        f"Recall after cores: {last_r_cores:.4f} | final recall: {last_r:.4f} | "
        f"true TEs: {len(true_te_set)} | "
        f"TP: {TP_list[-1] if TP_list else 0} | FP: {FP_list[-1] if FP_list else 0} | "
        f"TP_TE: {TP_te_list[-1] if TP_te_list else 0}"
    )

    norm = (
        LogNorm(vmin=c.min(), vmax=c.max())
        if (c > 0).all() and (c.max() / max(c.min(), 1e-12) > 10)
        else Normalize(vmin=c.min(), vmax=c.max())
    )
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    sc = ax.scatter(
        recall_list, precision_list, c=c, cmap="plasma", norm=norm, s=size, edgecolors="k"
    )
    ax.plot(recall_list, precision_list, color="grey", alpha=0.4, linewidth=1)
    ax.axvline(
        x=last_r_cores,
        color="green",
        linestyle="--",
        label="last recall with extended t-cores",
    )

    sizes_for_null = None
    if rate is not None:
        if core_sizes is not None:
            sizes_for_null = list(core_sizes) + [1] * n_unitigs
        else:
            sizes_for_null = [1] * n
        proba_comp = [size_to_proba(s, rate) for s in sizes_for_null]
        # Empty set: 100% precision, 0% recall
        recall_null = [0.0]
        precision_proba_cumul = [1.0]
        cumul = 0.0
        tot = 0
        for p in proba_comp:
            cumul += p
            tot += 1
            precision_proba_cumul.append(cumul / tot)
            recall_null.append(recall_list[tot - 1])
        ax.plot(
            recall_null,
            precision_proba_cumul,
            color="red",
            label="Cumulative theoretical probability",
        )

    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(title)
    ax.grid(True)
    ax.legend(loc="upper right")
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Max extended degree")
    plt.savefig(output_png, bbox_inches="tight")
    plt.show()
    plt.close()
    print(f"Saved: {output_png}")
    return precision_list, recall_list


def _read_true_te_file(path):
    true_te_set = set()
    with open(path, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            # allow TE_coverage_count_ab_filtered.txt (header + tab fields)
            if i == 0 and ("\t" in line or line.lower().startswith("te")):
                # skip header-like first line if it looks like a header
                parts = line.split("\t")
                if parts[0].lower() in {"te", "name", "id", "consensus"} or "count" in line.lower():
                    continue
            TE = line.split("\t")[0].strip()
            if TE:
                true_te_set.add(TE)
    return true_te_set


def _read_core_labels_file(path):
    labels = []
    with open(path, "r") as f:
        for line in f:
            # keep empty lines as false positives; strip newline only
            labels.append(line.rstrip("\n").rstrip("\r"))
    return labels


# ---------------------------------------------------------------------------
# Example list API (edit this block or import plot_pr_from_core_list)
# ---------------------------------------------------------------------------
EXAMPLE_CORE_TE_LIST = ["TE1","TE2", ""]#, "TE2", "TE1"]
EXAMPLE_TRUE_TE_SET = {"TE1", "TE2", "TE3"}
EXAMPLE_EXTENDED_DEGREE_CORES_LIST = [20,15,12]#,11,9]
EXAMPLE_SIZE_CORES_LIST = [5,10,6]#,7,3]#, 10, 3]

EXAMPLE_UNITIG_TE_LIST = [ "", "","", "TE3", ""]
EXAMPLE_EXTENDED_DEGREE_UNITIG_LIST = [9, 8, 6, 4, 2]

if __name__ == "__main__":
    Arg = sys.argv[:]

    # Minimal in-code example: python3 list_cores_roc_curve.py
    if len(Arg) == 1:
        if len(EXAMPLE_CORE_TE_LIST) != len(EXAMPLE_EXTENDED_DEGREE_CORES_LIST):
            raise ValueError(
                "EXAMPLE_CORE_TE_LIST and EXAMPLE_EXTENDED_DEGREE_CORES_LIST must have the same length"
            )
        if len(EXAMPLE_CORE_TE_LIST) != len(EXAMPLE_SIZE_CORES_LIST):
            raise ValueError(
                "EXAMPLE_CORE_TE_LIST and EXAMPLE_SIZE_CORES_LIST must have the same length"
            )
        if len(EXAMPLE_UNITIG_TE_LIST) != len(EXAMPLE_EXTENDED_DEGREE_UNITIG_LIST):
            raise ValueError(
                "EXAMPLE_UNITIG_TE_LIST and EXAMPLE_EXTENDED_DEGREE_UNITIG_LIST must have the same length"
            )
        plot_pr_from_core_list(
            EXAMPLE_CORE_TE_LIST,
            EXAMPLE_TRUE_TE_SET,
            output_png="pr_from_core_list_example.png",
            title="Precision-Recall Curve for TE Prediction",
            ranks=EXAMPLE_EXTENDED_DEGREE_CORES_LIST,
            core_sizes=EXAMPLE_SIZE_CORES_LIST,
            unitig_te_labels=EXAMPLE_UNITIG_TE_LIST,
            unitig_ranks=EXAMPLE_EXTENDED_DEGREE_UNITIG_LIST,
            rate=0.03,
        )
        exit(0)

    if len(Arg) not in [4, 5]:
        print(
            "Use :\n"
            f"  {Arg[0]}\n"
            f"  {Arg[0]} true_TEs.txt cores_labels.txt output.png [rate]\n"
            "\n"
            "Or import plot_pr_from_core_list([...], true_te_set, output_png=...)\n"
            'cores_labels.txt: one core per line, e.g. "TE1", "", "TE1; TE2"'
        )
        print("Arg used : ")
        print(Arg)
        exit(1)

    true_te_set = _read_true_te_file(Arg[1])
    core_te_labels = _read_core_labels_file(Arg[2])
    out = Arg[3]
    rate = float(Arg[4]) / 100.0 if len(Arg) == 5 else None

    print("Number of true TEs: " + str(len(true_te_set)))
    print("Number of cores: " + str(len(core_te_labels)))

    plot_pr_from_core_list(
        core_te_labels,
        true_te_set,
        output_png=out,
        rate=rate,
        core_sizes=[1] * len(core_te_labels) if rate is not None else None,
    )
