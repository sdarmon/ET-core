/*
 * Ce programme permet de loader un graphe de De Bruijn
 * et de calculer les poids des arêtes et des sommets à
 * un rayon donné en entrée.
 */

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include "graph.h"


// graph implementation
int main(int argc, char** argv)
{
    if (argc!=4 and argc!=6 and argc!=8 and argc!=10){
        cout << "Expected use of this program: \n\n\t" <<argv[0] << " file.nodes file.edges radius -k kmer -o output.txt -h ham_err\n" << endl;
        return 0;
    }

    vector<Edge> E;
    vector<Node> V;

    char* nodesPath = argv[1];
    char* edgesPath = argv[2];

    int ham_err = 0;

    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges,E);
    read_node_file(nodes,V);

    Graph G(V,E);
    if(argc >= 6 and argv[4][1]=='k' ){
        G.kmer = stoi(argv[5]);
    }
    if (argc >= 8 and argv[8][1]=='h' ){
        ham_err = stoi(argv[9]);
    }

    G.weighing();
    G.weighingAllNodes(stoi(argv[3]), ham_err);

    if(argc == 6 and argv[4][1]=='o'){
        ofstream output;
        output.open(argv[5]);
        printGraphVertices(G,output);
        output.close();
    } else if(argc == 8 and argv[6][1]=='o'){
        ofstream output;
        output.open(argv[7]);
        printGraphVertices(G,output);
        output.close();
    } else if(argc == 10 and argv[6][1]=='o'){
        ofstream output;
        output.open(argv[7]);
        printGraphVertices(G,output);
        output.close();
    } else {
        printGraphVertices(G);
    }
    return 0;

}