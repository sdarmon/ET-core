# Cette fonction permet de traiter les poids des sommets, après avoir
# utilisé l'executable graph.exe.
# L'option "dot" permet d'afficher les courbes
# L'option "his" permet d'afficher des histogrammes des données
# L'option "top10" permet d'obtenir le poids correspondant à la limite
# des 10% d'unitigs aux poids les plus élevés
# L'option "reverse" permet d'afficher les courbes avec les poids
# matchant avec un TE ou une répétition connu.

import numpy as np
import sys

Arg = sys.argv[:]

if len(Arg) not in [3, 4]:
    print("Wrong input used.\n Use : " + Arg[
        0] + " output.txt format optional_input.txt\n Where format is 'dot', 'his', 'top1', 'top10', 'top0001' or "
             "'reverse', the format of the plot.")

elif Arg[2] == "dot":

    import matplotlib.pyplot as plt

    y = []
    m = 0
    with open(Arg[1]) as f:
        for line in f:
            L = line.split('\t')
            if len(line) < 2 or int(L[2][:-1]) < 0:
                continue
            y.append(int(L[2][:-1]))
            m = max(m, y[-1])
        x = np.arange(1, m + 1, 1)
        Y = [0 for i in range(m)]
        for el in y:
            Y[el - 1] += 100 / len(y)
    xb = []
    yb = []

    xg = []
    yg = []

    for i in range(len(x)):
        if Y[i] == float(0):
            yb.append(Y[i])
            xb.append(x[i])
        else:
            yg.append(Y[i])
            xg.append(x[i])

    #fig, axs = plt.subplots(1)

    plt.plot(xg, yg, 'bo', label='Fréquence non nulle')
    plt.plot(x, Y, 'green')
    plt.ylabel("Frequences (en %)")
    plt.xlabel("Taille des poids")
    plt.title("Fréquence des différents poids (Pour un rayon de 10 nucléotides)")#fig.suptitle("Frequencies of the different sizes of weight (Case of " + Arg[0][:-4] + ")")


    #axs[1].plot(xb, yb, 'ro', label='Null Frequency')
    #axs[1].legend()
    #axs[1].plot(x, Y, 'green')
    #axs[1].set(ylabel="Frequences (en %)")


    plt.legend()
    plt.show()


elif Arg[2] == "his":
    import matplotlib.pyplot as plt

    y = []
    m = 0
    with open(Arg[1]) as f:
        for line in f:
            if len(line) < 2:
                continue
            poids = int(line.split("\t")[-1][:-1])
            m = max(poids,m)
            y.append(poids)
        x = np.arange(1, len(y) + 1, 1)

    Y = [0 for i in range(10)]
    X = [str(i) for i in range(10)]
    for i in range(len(y)):
        Y[(10 * y[i]) // (m+1)] +=1
    for i in range(10):
        Y[i] = 100 * Y[i] / len(y)
    plt.bar(X[1:],Y[1:], color='green')
    plt.ylabel("Fréquences (en %)")
    plt.xlabel("Taille des poids")
    plt.title("Histogramme des fréquences des différents poids par intervalle de 10% (premier intervalle non représenté)")#fig.suptitle("Frequencies of the different sizes of weight (Case of " + Arg[0][:-4] + ")")
    #plt.title("Histogram of the frequencies of the different sizes of weight (Case radius = 200)")
    plt.show()

elif Arg[2] == "top1":
    y = []
    m = 0
    size = 1
    with open(Arg[1]) as f:
        for line in f:
            L = line.split('\t')
            if len(line) < 2 or int(L[2][:-1]) < 0:
                continue
            y.append(int(L[2][:-1]))
            m = max(m, y[-1])
    x = np.arange(1, m + 1, 1)
    Y = [0 for i in range(m+1)]
    for el in y:
        Y[el] += 100 / len(y)
    M = 0
    for el in Y:
        M += el
        if M >= 99:
            break
        size += 1
    print(size)


elif Arg[2] == "top10":
    y = []
    m = 0
    size = 1
    with open(Arg[1]) as f:
        for line in f:
            L = line.split('\t')
            if len(line) < 2 or int(L[2][:-1]) < 0:
                continue
            y.append(int(L[2][:-1]))
            m = max(m, y[-1])
    x = np.arange(1, m + 1, 1)
    Y = [0 for i in range(m)]
    for el in y:
        Y[el - 1] += 100 / len(y)
    M = 0
    for el in Y:
        M += el
        if M >= 90:
            break
        size += 1
    print(size)

elif Arg[2] == "top0001":
    y = []
    m = 0
    size = 1
    with open(Arg[1]) as f:
        for line in f:
            L = line.split('\t')
            if len(line) < 2 or int(L[2][:-1]) < 0:
                continue
            y.append(int(L[2][:-1]))
            m = max(m, y[-1])
    x = np.arange(1, m + 1, 1)
    Y = [0 for i in range(m+1)]
    for el in y:
        Y[el] += 1000 / len(y)
    M = 0
    for el in Y:
        M += el
        if M >= 999:
            break
        size += 1
    print(size)

elif Arg[2] == "reverse":
    import matplotlib.pyplot as plt

    y = []
    m = 0
    ref = []
    Z = []
    with open(Arg[3]) as f:
        for line in f:
            L = line.split('\t')
            if len(line) < 2 or int(L[2][:-1]) < 0:
                continue
            ref.append(int(L[2][:-1]))
    with open(Arg[1]) as f:
        for line in f:
            if len(line) < 2:
                break
            L = line.split('\t')
            y.append(int(L[2][:-1]))
            m = max(m, y[-1])
    x = np.arange(1, m + 1, 1)
    Y = [0 for i in range(m)]
    Z = [0 for i in range(m)]
    for el in y:
        Y[el - 1] += 100 / len(y)
    for el in ref:
        Z[el - 1] += 100 / len(y)
    xb = []
    yb = []

    xg = []
    yg = []
    Z2 = [100 * el for el in Z]
    for i in range(len(x)):
        if Y[i] == float(0):
            yb.append(Y[i])
            xb.append(x[i])
        else:
            yg.append(Y[i])
            xg.append(x[i])

    fig, axs = plt.subplots(3)

    axs[0].plot(xg, yg, 'bo', label='Positive Frequency')
    axs[0].legend()
    axs[1].plot(xb, yb, 'ro', label='Null Frequency')
    axs[1].legend()
    axs[2].plot(x, Z, 'yo', label='Matching sequencies')
    axs[2].legend()
    axs[0].plot(x, Y, 'green')
    axs[1].plot(x, Y, 'green')
    axs[2].plot(x, Y, 'green')
    axs[2].plot(x, Z2, color='orange', linestyle='dashed', label="Frequencies times 100")
    axs[2].legend()
    axs[0].set(ylabel="Frequencies (in %)")
    axs[1].set(ylabel="Frequencies (in %)")
    axs[2].set(ylabel="Frequencies of matching sequencies (in %)")
    plt.xlabel("Sizes of weight")
    fig.suptitle("Frequencies of the different sizes of weight (Case radius = 200)")
    plt.show()

else:
    print("Wrong input used.\n Use : " + Arg[
        0] + "output.txt format optional_input.txt\n Where format is 'dot', 'his', 'top10' or 'reverse', the format "
             "of the plot.")
