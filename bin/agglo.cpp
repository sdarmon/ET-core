/*
 * Ce programme permet de faire une agglomeration des composantes
 * des graphes.
 */

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <set>
#include <iterator>
#include "graph.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstring-compare"
#define MAX 1024

// Permet de déterminer s'il existe un chemin de la composante i à la composante j, et s'il existe, renvoie la distance.
void chemin_local(int i, vector<int> &endings, Graph &G, vector<vector<int>> &components, int profondeurMax){
    vector<int> comp1;
    vector<Neighbor*> voisin1;
    vector<int> pronf;
    pronf.clear();
    voisin1.clear();
    comp1.clear();

    //Étape d'initialisation:
    for(int k= 0; k< components[i].size(); k++){
        comp1.push_back(components[i][k]); //On récupère les sommets de la composante que l'on souhaite étudier

        //On récupère les voisins de ces sommets
        for (vector<Neighbor>::iterator it = G.Neighbors(components[i][k])->begin(); it != G.Neighbors(components[i][k])->end(); ++it){
            if (find(comp1.begin(),comp1.end(), it->val) == comp1.end()){ 
                voisin1.push_back(&(*it));
                pronf.push_back(G.Vertices[it->val].label.size()-G.kmer+1);
            }
        }
    }
    //On regarde si la j-eme composante s'intersecte avec la comp1 (possible en fonction du critère de la composante
    //choisi.
    for (int j = 0; j<components.size(); j++){
        for(int k= 0; k< components[j].size(); k++){
            if (find(comp1.begin(),comp1.end(),components[j][k]) != comp1.end()){ 
                endings[j] = 0;
                continue;
            }
        }
    }

    int step = 1; //On effectue le BFS par couches successives : on traite à chaque boucle les sommets à distance `step`
    //de la composante
    int modif = 1; //Booléen indiquant s'il y a eu une modification à la `step`-ème étape du BFS
    int size;
    Neighbor* node;
    int pro; //Variable indiquant à quelle profondeur (en termes de nucléotides) de la composante

    while (modif && step < profondeurMax) { //Tant qu'un sommet a été rajouté à la couche précédente, on regarde tous les sommets de cette couche-là dans le BFS.
        modif = 0;
        size = voisin1.size();
        for (int k = 0; k<size; ++k){ //On boucle sur les sommets de la couche du BFS uniquement
            node = voisin1.front();
            pro = pronf.front();
            voisin1.erase(voisin1.begin());
            pronf.erase(pronf.begin());
            if (find(comp1.begin(),comp1.end(), node->val) != comp1.end()){ // Sommet déjà présent
                continue;
            }
            // Sommet également présent dans la composante d'arrivée
            for (int j = 0; j<components.size(); j++){
                if (find(components[j].begin(),components[j].end(),node->val) != components[j].end()){ 
                    if (endings[j] == -1){ //On garde toujours le plus court chemin
                        endings[j] = pro;
                    } else {
                        endings[j] = min(endings[j],pro);
                    }
                    continue;
                }
            }
            //Sinon on rajoute ce sommet comme vu et on ajoute ses voisins à traiter dans le BFS pour la prochaine couche
            comp1.push_back(node->val);
            modif = 1;
            for (vector<Neighbor>::iterator it = G.Neighbors(node->val)->begin(); it != G.Neighbors(node->val)->end(); ++it){
                if (it->label[0] != node->label[1]){
                    continue;
                }
                voisin1.push_back(&(*it));
                pronf.push_back(pro+G.Vertices[it->val].label.size()-G.kmer+1);
            }
        }
        step ++; //On passe à la couche suivante
    }
    //Aucun nouveau sommet n'a été rajouté à traiter... Il n'y a plus de chemins à visiter...
    return;
}

//Fonction permettant de récupérer l'index du sommet de `G` ayant le poids le plus élevé et qui n'a pas été déjà vu dans
//`vu_total` (i.e. `vu_total[index] == 0`)
int indexMax(Graph &G, vector<int> &vu_total){
    int index;
    int val = -1;
    for (int i = 0; i<G.N; i++){
        if (vu_total[i] == 0 && G.Vertices[i].weight > val){
            index = i;
            val = G.Vertices[i].weight;
        }
    }
    return index;
}

//Initialisation d'un vecteur avec la valeur -1. Il me semble que vector<int> vec(-1,n) doit fonctionner, à corriger!!
void initVec(vector<int> &vec, int n){
    vec.clear();
    for(int i = 0; i<n; i++){
        vec.push_back(-1);
    }
    return;
}

//Fonction permettant d'enregistrer une composante du graphe `G` dans un fichier
void save_comp(Graph &G, vector<int> &compo, string outputPrefix, int rang){
    ofstream output;
    output.open(outputPrefix+"/comp"+to_string(rang)+".txt");
    for (vector<int>::iterator it = compo.begin(); it != compo.end(); it++){
        output << *it << "\t" << G.Vertices[*it].label << "\t" << G.Vertices[*it].weight  << "\n";
    }
    return;
}

//Fonction indiquant si le setA est un sous-ensemble du setB. Le setA et le setB doivent être triés.
int subset(vector<int> &setA, vector<int> &setB){
    if (setA.size() > setB.size()){
        return 0;
    }
    int posi = 0;
    int posj = 0;
    while(posi< setA.size() and posj < setB.size()) {
        if (setA[posi] == setB[posj]){
            posi++;
            posj++;
        } else if (setA[posi] < setB[posj]){ //Cas où le i-ème élément n'a pas été trouvé dans le setB
            return 0;
        }
        else {
            posj++;
        }
    }
    return posi == setA.size();
}


int main(int argc, char** argv) {
    if (argc != 8) {
        cout << "Expected use of this program: \n\n\t" << argv[0]
             << " file.nodes file.edges -c value -k kmer outputPrefix \n" << endl;
        return 0;
    }

    //On charge le graphe
    vector <Edge> E;
    vector <Node> V;
    char *nodesPath = argv[1];
    char *edgesPath = argv[2];
    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges, E);
    read_node_file_weighted(nodes, V);

    Graph G(V, E);
    G.kmer = stoi(argv[6]);


    int index;
    //Construction du vector `vu_total`
    vector<int> vu_total(G.N, 0);

    //vector <vector<int>> components;
    //components.clear();

    //Ici le critère choisi pour définir une composante est juste d'être le sous-graphe connexe maximum du graphe G
    //auquel on a enlevé tous les sommets de poids inférieur au seuil `threshold`.
    int threshold;
    threshold = stoi(argv[4]);
    int weight;
    int compt;
    int m=0;
    vector<int> vu(G.N, 0);
    set<int> setVu;
    queue < Neighbor * > aVoir;

    index = indexMax(G, vu_total);
    weight = G.Vertices[index].weight;
    vu_total[index] = 1;

    while (weight >= threshold)
    {
        //Cleaning old variables
        for (int i = 0; i < G.N; i++) {
            vu[i] = 0;
        }
        vector<int> compo;
        compo.clear();
        while (!aVoir.empty()){
            aVoir.pop();
        }
        setVu.clear();

        //Initialisations
        vu[index] = 1;
        setVu.insert(index);
        for (vector<Neighbor>::iterator it = G.Neighbors(index)->begin(); it != G.Neighbors(index)->end(); ++it) {
            if (G.Vertices[it->val].weight >= threshold) {
                aVoir.push(&(*it));
            }
        }
        G.BFS_func(threshold, aVoir, vu,setVu);

        //Puis on ajoute les sommets trouvés à la composante
        compt=0;
        for (int i = 0; i < G.N; i++) {
            if (vu[i]) {
                compt++;
                compo.push_back(i);
                vu_total[i] = 1;
            }
        }
        if (compt > 1){
            save_comp(G, compo, argv[7], m); //On enregistre la composante sur l'ordinateur
            m++;
        }

        //Finalement, on recommence la boucle while
        index = indexMax(G, vu_total);
        weight = G.Vertices[index].weight;
        vu_total[index] = 1;

    }

    return 0;
}
