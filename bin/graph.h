#ifndef DEF_GRAPH
#define DEF_GRAPH

using namespace std;
#include <string>
#include <vector>
#include <set>
#include <queue>
#include <map>



// ===========================================================================
//                                  Structures
// ===========================================================================

/* 
 * Nœud d'un graphe.
 * Il est défini par une valeur val et par un label (string).
 * On rajoute aussi un argument weight pour lui donner un poids (variable)
 */
class Node {
public:
    int val, weight;
    string label;
    Node(int v, int w, string l);
};


/* 
 * Nœud dans la liste d'adjacence d'un graphe.
 * Il est défini par une valeur val et par un label (ici deux caractères correspondant au type de l'arête).
 * On rajoute aussi un argument weight pour lui donner un poids (variable)
 */
class Neighbor {
public:
    int val, weight;
    char label[3];
    Neighbor(int v, int w, char* l);
};

/* 
 * Liste des nœuds adjacents à un certain nœud donné.
 * Il est défini seulement par un vecteur de Neighbor.
 * Je pense qu'ici il serait préférable de ne pas utiliser une classe...
 */
class LstNode {
public:
    vector<Neighbor> adjVec;
    LstNode();
    LstNode(vector<Neighbor>& A);
};


/* 
 * Arête d'un graphe. 
 * Il est défini par une valeur de départ (start), d'arrivé (end) et un label (ici deux caractères correspondant au type de l'arête).
 * On rajoute aussi un argument weight pour lui donner un poids (variable)
 */
class Edge {
public:
    int start, end, weight;
    char label[3];
    Edge(int s, int e, int w, char* l);
};

typedef map<pair<int,int>,pair<int,int>> dic;


/* 
 * Arête d'un graphe. 
 * Il est défini par une valeur de départ (start), d'arrivé (end) et un label (ici deux caractères correspondant au type de l'arête).
 * On rajoute aussi un argument weight pour lui donner un poids (variable)
 */
class Graph
{
    public:

    //---------------
    // Variables
    //---------------

    int N;  // number of nodes in the graph
    int M; //number of edges in the graph
    int kmer = 41; //size of the kmers, by default = 41
    vector<LstNode> Adj;                //adjacency list as a vector of adjNodes
    vector<Node> Vertices;         //the vertices list as a vector of embNodes

    //---------------
    // Constructor
    //---------------

    Graph(vector<Node>& vertices, vector<Edge>& edges);

    //---------------
    // Methods
    //---------------
    bool operator() (int i,int j);

    void add(Edge &e);
    void add(int start,int end, int weight, char* labelEdge);
    void add(Node &v);
    void add(int val,int weight, string labelEdge);

    Neighbor* link(int start, int end);

    vector<Neighbor>* Neighbors(int n);
    void weighing();
    void BFS(int r, vector<Edge>& e ,vector<Neighbor*> &aVoir,vector<int> &vu);
    int BFSCount(vector<int> &rayons, int acc,vector<Neighbor*> &aVoir,vector<int> &vu);
    void BFS_func(int threshold, queue<Neighbor*> &aVoir,vector<int> &vu, set<int> & setVu);
    void BFS_func_dedupli(int threshold ,queue<int> &aVoir,vector<int> &vu, set<int> & setVu);
    void BFS_comp(vector<int> &seen,set<int> &vu, queue<Neighbor*> &aVoir, queue<int> &depth,queue<string> &labels,
                  vector<int> &sons, vector<int> &depthSons, vector<Neighbor*> &aretes, vector<string> & labelSons);
    void BFS_comp_GraphDedupli(vector<int> &seen,set<int> &vu, queue<Neighbor*> &aVoir, queue<int> &depth,queue<string> &labels,
                  vector<int> &sons, vector<int> &depthSons, vector<Neighbor*> &aretes, vector<string> & labelSons);
    void weighingANode(int source, int rayon);
    void weighingANodeGraphDuppli(int source, int rayon);
    void weighingAllNodes(int rayon, int ham_err);
    void weighingAllNodesGraphDuppli(int rayon);
    bool neigh(int u,int v,char s1,char s2);
    void BFScatch(vector<string> &kmers_at_distance_d,vector<Neighbor*> &aVoir,vector<int> &vu,vector<int> &rayons);
    int greedy_Hamming_cluster(vector<string> &kmers_at_distance_d, int d);
    void weighingANodeHamming(int source, int rayon, int d);
    };


//Other functions
void comp(char* s, char* r);
int reverse(int n, int N);
string reverse_complement(string a);
string make_label(string &a, string &b, char sensA, char sensB, int kmer);


void read_node_file( ifstream &node_file, vector<Node>& seqs);
void read_node_file_weighted( ifstream &node_file, vector<Node>& seqs);
void read_edge_file( ifstream &edge_file, vector<Edge>& edges );
void read_abundance_file( ifstream &ab_file, vector<double>& A );

void printGraphVertices(Graph& G);
void printGraphVertices(Graph& G,ofstream& output);
void printVertices(vector<Node>& V);
void printAbundance(vector<double> &A,ofstream& output);
void printVerticesBcalm(vector<Node>& V,ofstream& output);
void printGraphEdges(Graph& G);
void printEdges(vector<Edge>& E);
void printEdges(vector<Edge>& E,ofstream& output);
void printEdgesBcalm(vector<Edge>& E,ofstream& output);

#endif //def_graph