# Cette fonction permet de convertir les unitigs ayant un poids 
# supérieur au threshold au format fastq. Note à moi-même, cela 
# sera utilisé par STAR qui prend également en entrée les fichiers
# au format fasta qui sont donc bien moins gros et compliqué à faire!!!
# L'option "-reserve" permet de récupérer les unitigs qui sont issus du
# fichier d'intersectionKissNoDouble.txt


import sys

Arg = sys.argv[:]


def cleaning(sequence):
    """
    Entrée : une séquence de nucléotides
    Sortie : la séquence sans queue poly A de longueur au moins 15, en autorisant au plus une erreur
    """
    n = len(sequence)
    misstake = 1
    compt = 1
    last = sequence[0]
    while compt < n and (sequence[compt] == last or misstake):
        if sequence[compt] != last:
            misstake = 0
            if compt + 1 < n and sequence[compt + 1] != last:
                break
        compt += 1

    misstake = 1
    compt2 = 1
    last = sequence[n - 1]
    while compt2 < n and (sequence[n - compt2 - 1] == last or misstake):
        if sequence[n - compt2 - 1] != last:
            misstake = 0
            if n - compt2 - 2 < 0 and sequence[n - compt2 - 2] != last:
                break
        compt2 += 1

    if compt > 14 and compt2 > 14:
        return sequence[compt:n - compt2]
    elif compt > 14:
        return sequence[compt:]
    elif compt2 > 14:
        return sequence[:n - compt2]
    return sequence


def isPoly(sequence):
    """
    Entrée : une séquence de nucléotides
    Sortie : s'il existe une queue poly
    """
    if len(sequence) == 0:
        return True
    n = len(sequence)
    m = 0
    compteur = 1
    last = sequence[0]
    for i in range(1, n):
        if sequence[i] == last:
            compteur += 1
        else:
            m = max(m, compteur)
            compteur = 1
            last = sequence[i]
    m = max(m, compteur)
    if m > 14 or n == 0:
        return True
    return False


if len(Arg) not in [4, 5, 6]:
    print("Use : " + Arg[
        0] + "input.txt output.fq threshold [-reverse ref.txt] [-clean] \n Ne pas utiliser -reverse et -clean en même "
             "temps!")
    exit()
if len(Arg) == 4:
    seqs = []
    t = int(Arg[3])
    with open(Arg[1], 'r') as f:
        for line in f:
            if len(line) < 2:
                break
            L = line.split('\t')
            if '\n' in L[2]:
                th = L[2][:-1]
            else:
                th = L[2]
            if int(th) >= t:
                seqs.append(L[1])
    with open(Arg[2], 'w') as f:
        compt = 0
        for seq in seqs:
            if len(seq) == 0 :
                continue
            f.write(">SEQ_" + str(compt) + "\n")
            compt += 1
            f.write(seq + "\n")
elif len(Arg) == 5:  # i.e. argument -clean
    seqs = []
    t = int(Arg[3])
    with open(Arg[1], 'r') as f:
        for line in f:
            if len(line) < 2:
                break
            L = line.split('\t')
            if int(L[2][:-1]) >= t:
                seqs.append(line)
    with open(Arg[2], 'w') as f:
        compt = 0
        for seq in seqs:
            L = seq.split('\t')
            seqq = cleaning(L[1])
            if isPoly(seq):
                f.write(L[0] + "\t" + seqq + "\t-1\n")
            else:
                f.write(L[0] + "\t" + seqq + "\t" + L[2])
else: #i.e. argument -reverse
    seqs = []
    ref = set()
    t = int(Arg[3])
    with open(Arg[1], 'r') as f:
        for line in f:
            if len(line) < 2:
                break
            L = line.split('\t')
            if int(L[2][:-1]) >= t:
                seqs.append(line)
    with open(Arg[5], 'r') as f2:
        for line in f2:
            if len(line) < 2:
                break
            L = line.split('\t')
            ref.add(int(L[12].split('_')[1]))
    with open(Arg[2], 'w') as f:
        compt = 0
        for seq in seqs:
            seqq = cleaning(seq.split('\t')[1])
            if isPoly(seqq):
                continue
            if compt in ref:
                f.write(seq)
            compt += 1
