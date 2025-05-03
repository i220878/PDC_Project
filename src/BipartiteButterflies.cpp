#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

int main() {
    std::ifstream file("../input/fixed_graph.txt");
    // Line is first index of file
    std::string line;

    std::vector<int> offsetsV, edgesV;
    std::vector<int> offsetsU, edgesU;

    // Iterate over first 5 linse of the input, containing "AdjacencyHypergraph",
    // nv, mv, nu, mu
    // nv = Number of Vertices for bipartition V,
    // mv = Number of Edges for bipartition V,
    // nu = Number of Vertices for bipartition U,
    // mu = Number of Edges for bipartition U
    double nv, mv, nu, mu;
    getline(file, line);
    if (line != "AdjacencyHypergraph") {
        std::cerr << "Error: First line must be 'AdjacencyHypergraph'" << std::endl;
        return 1;
    }
    getline(file, line);
    std::istringstream(line) >> nv;
    getline(file, line);
    std::istringstream(line) >> mv;
    getline(file, line);
    std::istringstream(line) >> nu;
    getline(file, line);
    std::istringstream(line) >> mu;

    offsetsV.resize(nv);
    edgesV.resize(mv);
    
    offsetsU.resize(nu);
    edgesU.resize(mu);

    for (int i = 0; i < nv; ++i) {
        getline(file, line);
        std::istringstream(line) >> offsetsV[i];
    }
    for (int i = 0; i < mv; ++i) {
        getline()
    } 
}