// ===========================================================================
// ===========================================================================
//                               Classe de graphes
// ===========================================================================
// ===========================================================================


/* 
 * Ce fichier permet de manipuler les graphes pondérés avec des classes adaptées 
 * Pour voir les definitions des classes, allez dans le fichier graph.h.
 * Ce fichier contient seulement les fonctions et méthodes liées aux classes
 * définies dans le fichier graph.h.
 */

/* 
 * WARNING
 * All the data is coded with int data type. For graph with more than 2147483648
 * (> 2x10^9) nodes or edges, you will get trouble!
 * Instead of int, use unsigned int (for the same weight, extend this upper bound
 * to 4x10^9) or long long int (for twice the weight, extension to 2^63 ~ 10^19).
 *
 */

// ===========================================================================
//                               Include Libraries
// ===========================================================================
#include <list>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <sstream>
#include <queue>
// ===========================================================================
//                             Include Project Files
// ===========================================================================
#include "graph.h"


// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================
using namespace std;
#define MAX 1024


// ===========================================================================
//                                  Constructors
// ===========================================================================


Node::Node(int v, int w, string l){
        val = v;
        weight = w;
        label=(string)l;
    }

Neighbor::Neighbor(int v, int w, char* l){
        val = v;
        weight = w;
        strncpy(label, l, 2); // Copy max 2 characters
        label[2] = '\0';      // Ensure null-termination manually
    }

LstNode::LstNode() { // Je crois que ce constructeur ne fonctionne pas, à éviter
        adjVec.clear();
    }

LstNode::LstNode(vector<Neighbor>& A) {
        adjVec = A;
    }


Edge::Edge(int s, int e, int w, char* l) {
        start = s;
        end = e;
        weight = w;
        strcpy(label,l);
    }

/* Un graphe est défini par ses sommets (vecteur de Nodes) et
 * ses arêtes (des listes d'adjacence i.e. un vecteur de LstNode).
 * On tient également en mémoire le nombre de sommets (N) et le 
 * nombre d'arête M.
 */
Graph::Graph(vector<Node>& vertices, vector<Edge>& edges)  {
        N = vertices.size();
        M = 0;
        Vertices = vertices;
        Adj.clear();

        // initialize the LstNode for all vertices
        for (int i = 0; i < vertices.size(); i++){
            LstNode adj;
            Adj.push_back(adj);
            (Adj.back().adjVec).clear();
        }

        // construct directed graph by adding edges to it
        for (unsigned i = 0; i < edges.size(); i++)  {
            add(edges[i].start, edges[i].end, edges[i].weight, edges[i].label);
             }
    }

 bool Graph::operator() (int i,int j) { return (Vertices[i].weight <= Vertices[j].weight);}
//================================================================
//                  Public methods
//================================================================

//Ajoute une arête (de type Edge) au graphe.
void Graph::add( Edge &e) {
    Neighbor node(e.end,e.weight,e.label);
   (Adj[e.start].adjVec).push_back(node);
   M++;
    }

//Ajoute une arête (définie par ses quatre paramètres) au graphe.
void Graph::add(int start,int end, int weight, char* labelEdge)   {
        Neighbor node(end,weight,labelEdge);
       (Adj[start].adjVec).push_back(node);
       M++;
    }

//Ajoute un sommet (de type Node) au graphe.
void Graph::add(Node &v) {
   Vertices.push_back(v);
       LstNode A;
       Adj.push_back(A);
   N++;
    }

//Ajoute un sommet (défini par ses trois paramètres) au graphe.
void Graph::add(int val,int weight, string label)   {
        Node node(val,weight,label);
       Vertices.push_back(node);
       LstNode A;
       Adj.push_back(A);
       N++;
    }

/* Permet de savoir s'il existe une arête du sommet start au sommet end.
 * Si c'est le cas, revoit un pointeur vers l'arête dans la liste d'adjacence.
 * Sinon, revoit le pointeur NULL.
 */
Neighbor* Graph::link(int start, int end) {
    if (start >= N) {
        return NULL;
    }
    for (int index = 0; index < Adj[start].adjVec.size(); index++){
        if (Adj[start].adjVec[index].val == end) {
            return &Adj[start].adjVec[index];
        }
    }
    return NULL;
}

/* WARNING
 * Dans ce model de graphe, on suppose que l'on ne peut pas retirer des
 * sommets et donc la valeur d'un sommet est tout simplement son indice
 * dans le vecteur des sommets. La fonction suivante sert juste à
 * renvoyer la liste d'adjacence du n-ème sommet du graphe.
 * En cas d'ajout de la possibilité de retirer des sommets, cette
 * fonction sera bien à préciser car ne renverra pas forcément les
 * voisins d'un sommet de valeur donnée!!!!
 */ 
vector<Neighbor>* Graph::Neighbors(int n){
    return &(Adj[n].adjVec);
}


//Donne comme poids à toutes les arêtes la taille du label du sommet d'arrivée
void Graph::weighing(){
    for (int index = 0; index < N; index++){
        for (vector<Neighbor>::iterator it = Neighbors(index)->begin(); it != Neighbors(index)->end(); ++it){
            it->weight = Vertices[it->val].label.size();
        }
    }
}

//Calcule l'arête complémentaire de s.
void comp(char* s, char* r){
    if (s[0] == 'F')
    {
        r[1]='R';
    } else {
        r[1] = 'F';
    }
    if (s[1] == 'F')
    {
        r[0]='R';
    } else {
        r[0] = 'F';
    }
}

// BFS qui stocke les arêtes vues dans e. Le premier argument permet de dire dans quel sens on va
void Graph::BFS(int r, vector<Edge>& e ,vector<Neighbor*> &aVoir,vector<int> &vu){
    if (aVoir.size() == 0){ //Cas de terminaison, on a terminé le BFS
        return;
    }
    Neighbor* node = aVoir.front();
    aVoir.erase(aVoir.begin());
    if (find(vu.begin(),vu.end(),node->val) != vu.end()){ //Cas où le sommet a été vu par le BFS
        return BFS(r,e,aVoir,vu);
    }
    vu.push_back(node->val);
    for (vector<Neighbor>::iterator it = Neighbors(node->val)->begin(); it != Neighbors(node->val)->end(); ++it){
        //On boucle sur ses voisins
        if (it->label[0] == node->label[1]){
            //Cas où l'arrêt est bien valide, ce voisin est rajouté dans la file des visites
            aVoir.push_back(&(*it));
            if (r){
                Edge edge(node->val,it->val,0,it->label);
                e.push_back(edge);
            } else {
                char R[3];
                comp(it->label,R);
                Edge edge(it->val,node->val,0,R);
                e.push_back(edge);
            }
        }
    }
    return BFS(r,e,aVoir,vu); //Sinon, on continue sans prendre en compte le sommet.
}



// BFS qui teste si les sommets vérifient bien une condition (par une fonction)
void Graph::BFS_func(int threshold ,queue<Neighbor*> &aVoir,vector<int> &vu, set<int> & setVu){
    Neighbor* node ;
    while (aVoir.size() != 0){ //Cas de terminaison, on a terminé le BFS
        node = aVoir.front();
        aVoir.pop();
        if (vu[node->val]){ //Cas où le sommet a été vu par le BFS
            continue;
        }
        vu[node->val]=1;
        setVu.insert(node->val);

        for (vector<Neighbor>::iterator it = Neighbors(node->val)->begin(); it != Neighbors(node->val)->end(); ++it){
            //On boucle sur ses voisins
            if (Vertices[it->val].weight >= threshold and vu[it->val]==0){
                //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
                aVoir.push(&(*it));
            }
        }
    }
    return;
}
int reverse(int n, int N){
    if (n>= N/2){
        return n-N/2;
    } else {
        return n+N/2;
    }
}


// BFS qui teste si les sommets vérifient bien une condition (par une fonction)
void Graph::BFS_func_dedupli(int threshold ,queue<int> &aVoir,vector<int> &vu, set<int> & setVu){
    int vertex;
    while (aVoir.size() != 0){ //Cas de terminaison, on a terminé le BFS
        vertex = aVoir.front();
        aVoir.pop();
        if (vu[vertex]){ //Cas où le sommet a été vu par le BFS
            continue;
        }
        vu[vertex]=1;
        setVu.insert(vertex);

        for (vector<Neighbor>::iterator it = Neighbors(vertex)->begin(); it != Neighbors(vertex)->end(); ++it){
            //On boucle sur ses voisins
            if (Vertices[it->val].weight >= threshold and vu[it->val]==0){
                //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
                aVoir.push(it->val);
            }
        }
        for (vector<Neighbor>::iterator it = Neighbors(reverse(vertex,N))->begin(); it != Neighbors(reverse(vertex,N))->end(); ++it){
            //On boucle sur ses voisins
            if (Vertices[reverse(it->val,N)].weight >= threshold and vu[reverse(it->val,N)]==0){
                //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
                aVoir.push(reverse(it->val,N));
            }
        }
    }
    return;
}


typedef map<pair<int,int>,pair<int,int>> dic;

string reverse_complement(string a){
    string aux = "";
    for (int i = a.size()-1; i>=0; i--){
        if (a[i] == 'A'){
            aux = aux + 'T';
        } else if (a[i] == 'T'){
            aux = aux + 'A';
        } else if (a[i] == 'G'){
            aux = aux + 'C';
        } else {
            aux = aux + 'G';
        }
    }
    return aux;
}

int min(int a, int b){
    if (a < b) {
        return a;
    }
    return b;
}

string make_label(string &a, string &b, char sensA, char sensB, int kmer){ //Voir explication dans le cahier
    string aux;
    if( sensA == 'F' and sensB == 'F'){
        aux = a + b.substr(min(kmer-1,b.size()-1));
    } else if (sensA == 'R' and sensB == 'R'){
        aux = a + reverse_complement(b).substr(min(kmer-1,b.size()-1));
    } else if (sensA == 'F' and sensB == 'R'){
        aux = a + reverse_complement(b).substr(min(kmer-1,b.size()-1));
    } else {
        aux = a + b.substr(min(kmer-1,b.size()-1));
    }
    return aux;
}

void Graph::BFS_comp(vector<int> &seen,set<int> &vu, queue<Neighbor*> &aVoir, queue<int> &depth, queue<string> &labels,
                     vector<int> &sons, vector<int> &depthSons, vector<Neighbor*> &aretes, vector<string> & labelSons){
    int pronf;
    string label;
    while (aVoir.size() != 0) { //Cas de terminaison, on a terminé le BFS
        Neighbor *node = aVoir.front();
        aVoir.pop();
        pronf = depth.front();
        depth.pop();
        label = labels.front();
        labels.pop();
        if (vu.find(node->val) != vu.end()) { //Cas où le sommet a été vu par le BFS
            //Cela veut dire qu'on l'a vu par un autre chemin, et il y a nécessairement une bulle par ici
            continue;
        }
        vu.insert(node->val);
        if (seen[node->val] <= 0) { //Cas où le sommet est en dehors du périmètre, on sort
            if(Vertices[node->val].label.size()>=kmer){
                sons.push_back(node->val);
                aretes.push_back(node);
                depthSons.push_back(pronf);
                labelSons.push_back(label);
            }
            continue;
        }
        if (Vertices[node->val].label.size() < kmer) {
            continue;
        }
        for (vector<Neighbor>::iterator it = Neighbors(node->val)->begin(); it != Neighbors(node->val)->end(); ++it) {
            //On boucle sur ses voisins
            if (it->label[0] == node->label[1]) {
                //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
                aVoir.push(&(*it));
                depth.push(pronf+1);
                labels.push(make_label(label,Vertices[node->val].label,node->label[0],node->label[1],kmer));
            }
        }
    }
    return;
}


void Graph::BFS_comp_GraphDedupli(vector<int> &seen,set<int> &vu, queue<Neighbor*> &aVoir, queue<int> &depth, queue<string> &labels,
                                  vector<int> &sons, vector<int> &depthSons, vector<Neighbor*> &aretes, vector<string> & labelSons){
    int pronf;
    string label;
    while (aVoir.size() != 0) { //Cas de terminaison, on a terminé le BFS
        Neighbor *node = aVoir.front();
        aVoir.pop();
        pronf = depth.front();
        depth.pop();
        label = labels.front();
        labels.pop();
        if (vu.find(node->val) != vu.end()) { //Cas où le sommet a été vu par le BFS
            //Cela veut dire qu'on l'a vu par un autre chemin, et il y a nécessairement une bulle par ici
            continue;
        }
        vu.insert(node->val);
        if (seen[node->val] <= 0) { //Cas où le sommet est en dehors du périmètre, on sort
            if(Vertices[node->val].label.size()>=kmer){
                sons.push_back(node->val);
                aretes.push_back(node);
                depthSons.push_back(pronf);
                labelSons.push_back(label);
            }
            continue;
        }
        if (Vertices[node->val].label.size() < kmer) {
            continue;
        }
        for (vector<Neighbor>::iterator it = Neighbors(node->val)->begin(); it != Neighbors(node->val)->end(); ++it) {
            //On boucle sur ses voisins
            if (it->label[0] == node->label[1]) {
                //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
                aVoir.push(&(*it));
                depth.push(pronf+1);
                labels.push(label+Vertices[node->val].label.substr(kmer-1));
            }
        }
    }
    return;
}

/* /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ 
 *   Petite approximation ici, je suppose que dès lors que l'on 
 *   voit un sommet par le sens forward, on ne le recroisera pas
 *   par le sens reverse ! Cela semble peu probable d'arriver et
 *   l'impact semble peu important (car BFS) mais on sous-estime
 *   la vraie valeur.
 *  /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ 
 *
 * Fonction inductive permettant d'effectuer un BFS, en mettant à jour
 * les données pour les appels inductifs suivants. De manière générale,
 * cette fonction sert à compter le nombre de sommets dans une boule de
 * rayon fixé, autour d'un sommet donné. 
 */
int Graph::BFSCount(vector<int> &rayons, int acc,vector<Neighbor*> &aVoir,vector<int> &vu){
    if (aVoir.size() == 0){ //Cas de terminaison, on a terminé le BFS
        return acc;
    }
    int rayon = rayons.front();
    rayons.erase(rayons.begin());
    Neighbor* node = aVoir.front();
    aVoir.erase(aVoir.begin());
    if (find(vu.begin(),vu.end(),node->val) != vu.end()){ //Cas où le sommet a été vu par le BFS
        return BFSCount(rayons,acc,aVoir,vu);
    }
    vu.push_back(node->val);
    if (node->weight -kmer+1 <= rayon) { //Cas où le sommet est bien dans la boule
        int nb_neigh = Neighbors(node->val)->size();
        for (vector<Neighbor>::iterator it = Neighbors(node->val)->begin(); it != Neighbors(node->val)->end(); ++it) {
            //On boucle sur ses voisins
            if (it->label[0] == node->label[1] && find(vu.begin(), vu.end(), it->val) == vu.end()) {
                //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
                aVoir.push_back(&(*it));
                rayons.push_back(rayon + kmer - 1 - node->weight);
            }
        }
        if (nb_neigh == 2) {
            return BFSCount(rayons, acc, aVoir, vu); //Cas d'un chemin, on prend en compte le sommet
        } else {
            return BFSCount(rayons, acc + 1, aVoir, vu); //On traite les cas suivants, en prenant en compte le sommet
        }
    }

    return BFSCount(rayons,acc,aVoir,vu); //Sinon, on continue sans prendre en compte le sommet.
}

//Permet de donner un poids à un sommet, correspondant aux nombres de sommets présents dans un rayon donné.
void Graph::weighingANode(int source, int rayon) {
    //Initialisation des files pour le bfs. Il serait plus judicieux ici d'utiliser un set<int> pour vu!!!!
    vector<Neighbor*> aVoir;
    vector<int> vu;
    vector<int> rayons; //Ici on a une file de rayons car comme on calcule la distance en termes de nucléotides, cette
                        // métrique est différente de la distance entre sommets dans le graphe compacté

    //De plus, pour la définition de la boule autour d'un unitig, il faut prendre en compte le fait que ce dernier est
    //peut-être bien plus grand qu'un k-mer! Le critère choisi pour prendre en compte ce paramètre est de pondérer les
    //rayons en fonction de toutes les positions possibles dans le graphe de De Bruijn non compacté et de choisir la
    //valeur de poids minimal pour bloquer gros unitigs. Voir cahier.
    int mini = 0;
    int firstTime = 1;
    int taille = Vertices[source].label.size() - kmer+1;
    for(int position = 0; position < taille; position++){ //On boucle sur toutes les positions
        rayons.clear();
        aVoir.clear();
        vu.clear();
        vu.push_back(source);
        for (vector<Neighbor>::iterator it = Neighbors(source)->begin(); it != Neighbors(source)->end(); ++it){
            if (it->label[0] == 'F'){
                aVoir.push_back(&(*it));
                rayons.push_back(rayon-position);
            } else {
                aVoir.push_back(&(*it));
                rayons.push_back(rayon-taille+1+position);
            }
        }
        if (firstTime) {
            mini = BFSCount(rayons,1,aVoir,vu);
            firstTime = 0;
        } else {
            mini = min(mini,BFSCount(rayons,1,aVoir,vu)); // On garde que le minimum des valeurs trouvées
        }
    }
    Vertices[source].weight = mini;
}

//Definition de BFScatch permettant de récupérer les strings des sommets à distance le rayon du sommet source
void Graph::BFScatch(vector<string> &kmers_at_distance_d,vector<Neighbor*> &aVoir,vector<int> &vu,vector<int> &rayons) {
    if (aVoir.size() == 0) { //Cas de terminaison, on a terminé le BFS
        return;
    }
    int rayonCourant = rayons.front();
    rayons.erase(rayons.begin());
    Neighbor* node = aVoir.front();
    aVoir.erase(aVoir.begin());


    if (rayonCourant == 0) {
        if (find(vu.begin(),vu.end(),node->val) != vu.end()) { //Cas où le sommet a été vu par le BFS
            vu.push_back(node->val);
            if (node->label[1] == 'F') {
                //On ajoute le 1er k-mer de node
                kmers_at_distance_d.push_back(Vertices[node->val].label.substr(0,kmer));
            } else {
                //On ajoute le 1er k-mer du reverse complement de node
                kmers_at_distance_d.push_back(reverse_complement(Vertices[node->val].label).substr(0,kmer));
            }
            return BFScatch(kmers_at_distance_d,aVoir,vu,rayons);
        }
    }

    if (find(vu.begin(),vu.end(),node->val) != vu.end()) { //Cas où le sommet a été vu par le BFS
        return BFScatch(kmers_at_distance_d,aVoir,vu,rayons);
    }
    vu.push_back(node->val);
    if (rayonCourant + kmer - 1 >= Vertices[node->val].label.size()) {
        for (vector<Neighbor>::iterator it = Neighbors(node->val)->begin(); it != Neighbors(node->val)->end(); ++it){
            if (it->label[0] == node->label[1]){
                aVoir.push_back(&(*it));
                rayons.push_back(rayonCourant + kmer - 1 - Vertices[node->val].label.size());
            }
        }
    } else {
        //On ajoute le k-mer de node
        if (node->label[1] == 'F') {
            kmers_at_distance_d.push_back(Vertices[node->val].label.substr(rayonCourant,kmer));
        } else {
            kmers_at_distance_d.push_back(reverse_complement(Vertices[node->val].label).substr(rayonCourant,kmer));
        }
    }
    return BFScatch(kmers_at_distance_d,aVoir,vu,rayons);
}

//Definition de greedy_Hamming_cluster permettant de trouver le nombre de classes de sequences de kmers_at_distance_d
//qui sont à distance d'Hamming d'au plus d nucléotides de l'une des autres séquences de la classe.
int Graph::greedy_Hamming_cluster(vector<string> &kmers_at_distance_d, int d){
    int taille = kmers_at_distance_d.size();
    int nb_classes = 0;
    vector<int> classes(taille,0);
    for (int i = 0; i < taille; i++){
        if (classes[i] == 0){
            nb_classes++;
            classes[i] = nb_classes;
        }
        for (int j = i+1; j < taille; j++){
            if (classes[j] == 0){
                int distance = 0;
                for (int k = 0; k < kmers_at_distance_d[i].size(); k++){
                    if (kmers_at_distance_d[i][k] != kmers_at_distance_d[j][k]){
                        distance++;
                    }
                }
                if (distance <= d){
                    classes[j] = classes[i];
                }
            } else {
                int distance = 0;
                for (int k = 0; k < kmers_at_distance_d[i].size(); k++){
                    if (kmers_at_distance_d[i][k] != kmers_at_distance_d[j][k]){
                        distance++;
                    }
                }
                if (distance <= d){
                    int ci = classes[i];
                    int cj = classes[j];
                    int cm = min(ci,cj);
                    for (int k = 0; k < taille; k++){
                        if (classes[k] == ci or classes[k] == cj){
                            classes[k] = cm;
                        }
                    }
                }
            }
        }
    }
    return nb_classes;
}

//Permet de donner un poids à un sommet, correspondant aux nombres de sommets à distance le rayon qui sont distants (distance de Hamming) de moins de d nucléotides.
void Graph::weighingANodeHamming(int source, int rayon, int d){
    vector<string> kmers_at_distance_d;
    kmers_at_distance_d.clear();
    vector<Neighbor*> aVoir;
    vector<int> vu;
    vector<int> rayons;
    int taille_hamming  ;
    //On rajoute dans Avoir les voisins du sommet source
    vu.push_back(source);
    for (vector<Neighbor>::iterator it = Neighbors(source)->begin(); it != Neighbors(source)->end(); ++it){
        if (it->label[0] == 'F'){
            aVoir.push_back(&(*it));
            rayons.push_back(rayon-1);
        } else {
            aVoir.push_back(&(*it));
            rayons.push_back(rayon-1);
        }
    }
    //On récupère les strings des sommets à distance rayon du sommet source
    BFScatch(kmers_at_distance_d,aVoir,vu,rayons);

    //On calcule le nombre de classes de séquences de kmers_at_distance_d qui sont à distance d'Hamming d'au plus d nucléotides de l'une des autres séquences de la classe.
    taille_hamming = greedy_Hamming_cluster(kmers_at_distance_d,d);

    //On donne ce nombre de classes comme poids au sommet source
    Vertices[source].weight = taille_hamming;
}


//Permet de donner un poids à un sommet, correspondant aux nombres de sommets présents dans un rayon donné.
void Graph::weighingANodeGraphDuppli(int source, int rayon) {
    //Initialisation des files pour le bfs. Il serait plus judicieux ici d'utiliser un set<int> pour vu!!!!
    vector<Neighbor*> aVoir;
    vector<int> vu;
    vector<int> rayons; //Ici on a une file de rayons car comme on calcule la distance en termes de nucléotides, cette
    // métrique est différente de la distance entre sommets dans le graphe compacté
    rayons.clear();
    aVoir.clear();
    vu.clear();
    vu.push_back(source);
    for (vector<Neighbor>::iterator it = Neighbors(source)->begin(); it != Neighbors(source)->end(); ++it){
        aVoir.push_back(&(*it));
        rayons.push_back(rayon-(Vertices[source].label.size() - kmer));
    }
    Vertices[source].weight = BFSCount(rayons,1,aVoir,vu);
}

//Permet de donner un poids à tous les sommets du graphe.
void Graph::weighingAllNodes(int rayon) {
    for (vector<Node>::iterator it = Vertices.begin(); it != Vertices.end(); ++it){
        weighingANode(it->val, rayon);
    }
}

//Permet de donner un poids à tous les sommets du graphe.
void Graph::weighingAllNodesGraphDuppli(int rayon) {
    for (vector<Node>::iterator it = Vertices.begin(); it != Vertices.end(); ++it){
        weighingANodeGraphDuppli(it->val, rayon);
    }
}

bool Graph::neigh(int u,int v,char s1,char s2){
    for (vector<Neighbor>::iterator it = Neighbors(u)->begin(); it != Neighbors(u)->end(); ++it) {
        if (it->val == v and it->label[0] == s1 and it->label[1] == s2) {
            return true;
        }
    }
    return false;
}

//================================================================
//                  Reading functions
//================================================================


//Lit un fichier contenant les sommets du graphe et les ajoute au vecteur seqs
void read_node_file( ifstream &node_file, vector<Node>& seqs)
{
    seqs.clear();
    string line;
    int u;
    string p;
    int compt;
    while (getline(node_file, line)) {
        istringstream ss(line);
        string substr;
        compt=0;
        while (getline(ss, substr, '\t')) {
            if (compt == 0) {
                u = stoi(substr);
            }
            else if (compt == 1){
                p=substr;
            }
            compt++;
        }
        Node node(u,0,p);
        seqs.push_back(node);
    }
}

//Lit un fichier contenant les sommets du graphe et les ajoute au vecteur seqs
void read_node_file_weighted( ifstream &node_file,vector<Node>& seqs)
{
    seqs.clear();
    string line;
    int u,v;
    string p;
    int compt;
    while (getline(node_file, line)) {
        istringstream ss(line);
        string substr;
        compt=0;
        while (getline(ss, substr, '\t')) {
            if (compt == 0) {
                u = stoi(substr);
            }
            else if (compt == 1){
                p=substr;
            } else if (compt == 2){
                v = stoi(substr);
            }
            compt++;
        }
        Node node(u,v,p);
        seqs.push_back(node);
    }
    return;
}


//Lit un fichier d'arêtes et les ajoute au vecteur edges
void read_edge_file( ifstream &edge_file, vector<Edge>& edges ) {
    edges.clear();
    string line;
    int u,v;
    char p[3];
    int compt;
    while (getline(edge_file, line)) {
        istringstream ss(line);
        string substr;
        compt=0;
        while (getline(ss, substr, '\t')) {
            if (compt == 0) {
                u = stoi(substr);
            }
            else if (compt == 1){
                v = stoi(substr);
            } else if (compt == 2){
                strcpy( p, substr.c_str() );
            }
            compt++;
        }
        Edge e(u,v,0,p);
        edges.push_back(e);
    }
}

void read_abundance_file( ifstream &ab, vector<double>& A ) {
    string line;
    while (getline(ab, line)) {
        A.push_back(stod(line));
    }
    return;
}

//================================================================
//                  Fonctions d'affichage
//================================================================



//Affiche les sommets d'un graphe
void printGraphVertices(Graph& G)
{
    vector<Node> V = G.Vertices;
    for (vector<Node>::iterator it = V.begin(); it != V.end(); ++it) {
        cout << it->val << "\t" << (string)it->label << "\t" << it->weight << "\n";
    }
    cout << endl;
}

//Affiche les sommets d'un graphe dans une autre sortie (dans un fichier par exemple)
void printGraphVertices(Graph& G,ofstream& output)
{
    vector<Node> V = G.Vertices;
    for (vector<Node>::iterator it = V.begin(); it != V.end(); ++it) {
        output << it->val << "\t" <<(string)it->label << "\t" << it->weight << "\n";
    }
}

//Affiche l'abondance dans une autre sortie (dans un fichier par exemple)
void printAbundance(vector<double> &A,ofstream& output)
{
    for (int i = 0; i < A.size(); i++) {
        output <<A[i] << "\n";
    }
    return;
}

//Affiche les sommets d'un graphe à partir d'un vecteur de Nodes
void printVertices(vector<Node>& V)
{
    for (vector<Node>::iterator it = V.begin(); it != V.end(); ++it) {
        cout << it->val << "\t" << (string)it->label << "\t" << it->weight << "\n";
    }
    cout << endl;
}

//Affiche les sommets d'un graphe à partir d'un vecteur de Nodes sans les poids
void printVerticesBcalm(vector<Node>& V,ofstream& output)
{
    for (vector<Node>::iterator it = V.begin(); it != V.end(); ++it) {
        output << it->val << "\t" << (string)it->label << "\n";
    }
    cout << endl;
}

//Affiche les arêtes d'un graphe
void printGraphEdges(Graph& G)
{
    vector<LstNode> Adj = G.Adj;
    int N = G.N;
    for (int i = 0; i<N ; i++ ) {
        for (vector<Neighbor>::iterator it = Adj[i].adjVec.begin(); it != Adj[i].adjVec.end(); ++it) {
            cout << i << " to " << it->val << " of type " << (string)it->label << " and of weight " << it->weight << endl;
        }
    }
    cout << endl;
}

//Affiche les arêtes d'un graphe à partir d'un vecteur d'Edges
void printEdges(vector<Edge>& E)
{
    int M = E.size();
    Edge* e;
    for (int i = 0; i<M ; i++ ) {
        *e = E[i];
        cout << e->start << " to " << e->end << " of type " << (string)e->label << " and of weight " << e->weight << endl;
    }
    cout << endl;
}

//Affiche les arêtes d'un graphe à partir d'un vecteur d'Edges
void printEdges(vector<Edge>& E,ofstream& output)
{
    for (vector<Edge>::iterator it = E.begin(); it != E.end(); ++it) {
        output << it->start << "\t" << it->end << "\t" << (string)it->label << "\t" << it->weight << "\n";
    }
}

//Affiche les arêtes d'un graphe à partir d'un vecteur d'Edges sans le poids
void printEdgesBcalm(vector<Edge>& E,ofstream& output)
{
    for (vector<Edge>::iterator it = E.begin(); it != E.end(); ++it) {
        output << it->start << "\t" << it->end << "\t" << (string)it->label << "\n";
    }
}