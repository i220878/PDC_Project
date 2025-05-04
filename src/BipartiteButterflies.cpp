#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <numeric>
using namespace std;

template<typename T>
void printOffsetsAndEdges(std::vector<T> off1,
std::vector<T> edge1, std::vector<T> off2, std::vector<T> edge2) {
    std::cout << "Offsets for bipartition V:\n";
    for (const auto& offset : off1) {
        std::cout << offset << ' ';
    }
    std::cout << "\nEdges for bipartition V:\n";
    for (const auto& edge : edge1) {
        std::cout << edge << ' ';
    }
    std::cout << "\nOffsets for bipartition U:\n";
    for (const auto& offset : off2) {
        std::cout << offset << ' ';
    }
    std::cout << "\nEdges for bipartition U:\n";
    for (const auto& edge : edge2) {
        std::cout << edge << ' ';
    }
    std::cout << '\n';
}

void printAdjacencyList(int** adj, int v1, int v2) {
    int n = v1 + v2;
    // width to print each index (for alignment)
    int w = 2; // adjust if your n > 99

    // 1) Column headers
    //    indent to align under row-labels (w+1 chars)
    bool recount = false;
    cout << string(w+3, ' ');
    for (int j = 0; j < n; ++j) {
        if (recount) cout << setw(w) << j - v1 << ' ';
        else cout << setw(w) << j << ' ';
        if (j == v1-1){
            cout << "| ";
            recount = true;
        }
    }
    cout << "\n";
    // 2) Separator line
    cout << string(w+1, ' ');
    for (int j = 0; j < n; ++j) {
        cout << string(w, '-') << '-';
        if (j == v1-1) cout << "-";
    }
    cout << "\n";

    // 3) Each row
    recount = false;
    for (int i = 0; i < n; ++i) {
        // row label
        if (recount) cout << setw(w) << i - v1 << ' ';
        else cout << setw(w) << i << ' ';
        cout << "| ";
        if (i == v1 - 1) {
            recount = true;
        }

        // row entries
        for (int j = 0; j < n; ++j) {
            cout << setw(w) << adj[i][j] << ' ';
            if (j == v1-1) cout << "| ";
        }
        cout << "\n";
    }
}

int main() {
    std::ifstream file("../input/fixed_graph.txt");
    // Line is first index of file
    std::string line;

    std::vector<int> offsetsV, edgesV;
    std::vector<int> offsetsU, edgesU;

    // Iterate over first 5 linse of the input, containing "AdjacencyHypergraph",
    // nv, mv, nu, mu
    // vV = Number of Vertices for bipartition V,
    // eV = Number of Edges for bipartition V,
    // vU = Number of Vertices for bipartition U,
    // eU = Number of Edges for bipartition U
    int v1, e1, v2, e2;
    getline(file, line);
    if (line != "AdjacencyHypergraph") {
        std::cerr << "Error: First line must be 'AdjacencyHypergraph'" << std::endl;
        return 1;
    }
    getline(file, line);
    std::istringstream(line) >> v1;
    getline(file, line);
    std::istringstream(line) >> e1;
    getline(file, line);
    std::istringstream(line) >> v2;
    getline(file, line);
    std::istringstream(line) >> e2;

    offsetsV.resize(v1);
    edgesV.resize(e1);
    
    offsetsU.resize(v2);
    edgesU.resize(e2);

    for (int i = 0; i < v1; ++i) {
        getline(file, line);
        std::istringstream(line) >> offsetsV[i];
    }
    for (int i = 0; i < e1; ++i) {
        getline(file, line);
        std::istringstream(line) >> edgesV[i];
    } 
    for (int i = 0; i < v2; ++i) {
        getline(file, line);
        std::istringstream(line) >> offsetsU[i];
    }
    for (int i = 0; i < e2; ++i) {
        getline(file, line);
        std::istringstream(line) >> edgesU[i];
    }

    // Vertices for V + Vertices for U for the adjacency list
    int** adj = new int*[v2 + v1];
    for (int i = 0; i < v2 + v1; ++i) {
        adj[i] = new int[v2 + v1];
        for (int j = 0; j < v2 + v1; ++j)
            adj[i][j] = 0;
    }
    
    // 0..vV - 1 rows and cols are both for first bipartition
    // vV..vU - 1 rows and cols are both for second bipartition
    for (int i = 0; i < v1; ++i) {
        int startV = offsetsV[i];
        // since the end offset is non-inclusive,
        // and if it's the last vertex, then it has as many edges
        // as total left to occupy
        int endV = i + 1 < v1 ? offsetsV[i + 1] : e1;
        
        for (int j = startV; j < endV; ++j) {
            int edge = edgesV[j];
            adj[i][v1 + edge] = 1;
            adj[v1 + edge][i] = 1;
        }
    }

    for (int i = 0; i < v2; ++i) {
        int startV = offsetsU[i];
        int endV = i + 1 < v2 ? offsetsU[i + 1] : e2;
        
        for (int j = startV; j < endV; ++j) {
            int edge = edgesU[j];
            adj[i + v1][edge] = 1;
            adj[edge][i + v1] = 1;
        }
    }
    

    /*
        STEP 1: COMPUTE VERTEX RANK -- USING DEGREE ORDER
    */
    int totalVertices = v1 + v2;
    vector<int> degrees(totalVertices);

    for (int i = 0; i < totalVertices; ++i) {
        int curDegree = 0;
        for (int j = 0; j < totalVertices; ++j) {
            if (adj[i][j] == 1) curDegree += 1;
        }
        degrees[i] = curDegree;
    }

    vector<int> sortedVertices(totalVertices);
    for (int i = 0; i < totalVertices; ++i) {
        sortedVertices[i] = i;
    }

    // Sort vertices in increasing (ascending) order of rank
    std::sort(sortedVertices.begin(), sortedVertices.end(), [&degrees](int a, int b) {
        return degrees[a] < degrees[b];
    });

    vector<int> ranks(totalVertices);

    for (int i = 0; i < totalVertices; ++i) {
        int idx = sortedVertices[i];
        ranks[idx] = i;
    }

    vector<vector<int>> neighbors(totalVertices);

    for (int i = 0; i < totalVertices; ++i) {
        for (int j = 0; j < totalVertices; ++j) {
            if (adj[i][j] == 1)
                neighbors[i].push_back(j);
        }
    }

    // Sorting neighbours in decreasing (descending) order of rank[neighbour]
    for (int i = 0; i < totalVertices; ++i) {
        std::sort(neighbors[i].begin(), neighbors[i].end(), [&ranks](int a, int b) {
            return ranks[a] > ranks[b];
        });
    }

    vector<vector<int>> modifiedNeighbors(totalVertices);

    /*
        We also define a modified neighborhood,
        Nv (u), to be the set of neighbors u′ ∈ N (u) such that
        rank(u′) > rank(v)
    */
    for (int v1 = 0; v1 < totalVertices; ++v1) {
        for (int v2 : neighbors[v1]) {
            if (ranks[v2] > ranks[v1]) {
                modifiedNeighbors[v1].push_back(v2);
            }
        }
    }
    /*
        We define a modified degree,
        degv (u), to be the number of neighbors u′ ∈ N (u) such that
        rank(u′) > rank(v

        NOTE: MODIFIED DEGREE is simply modifiedNeighbors[index].size()
    */

    /*
        STEP 2: WEDGE RETRIEVAL / ENUMERATION
    */// Build modifiedNeighbors based on ranks
    // 2: Use PREFIX-SUM to compute a function I that maps wedges to indices in order
    // Build immediateNeighbors: exclusive prefix‐sum of modified‐degrees
    vector<int> immediateNeighbors(totalVertices + 1);
    immediateNeighbors[0] = 0;
    for (int i = 0; i < totalVertices; ++i) {
        immediateNeighbors[i + 1] =
            immediateNeighbors[i] + (int)modifiedNeighbors[i].size();
    }
    int M = immediateNeighbors[totalVertices];

    // 2) For each (u1, v), we calculate how many possible u2's it results in, i.e how many neighbors of v that have a rank
    // greater than v

    vector<int> secondNeighbors(M);
    for (int u1 = 0; u1 < totalVertices; ++u1) {
        int modDeg = (int)modifiedNeighbors[u1].size();
        for (int i = 0; i < modDeg; ++i) {
            int v = modifiedNeighbors[u1][i];
            // count *all* neighbors of v except the back‐edge to u1
            int cnt2 = 0;
            for (int w : neighbors[v]) {
                if (ranks[w] > ranks[u1]) {
                    ++cnt2;
                }
            }
            secondNeighbors[immediateNeighbors[u1] + i] = cnt2;
        }
    }

    // 3) Calculate exclusive prefix‐sum of 2nd modified degrees (neighbors of neighbors of u in )
    vector<size_t> secondNeighborsPrefixSum(M + 1);
    secondNeighborsPrefixSum[0] = 0;
    for (int k = 0; k < M; ++k) {
        secondNeighborsPrefixSum[k + 1] =
            secondNeighborsPrefixSum[k] + secondNeighbors[k];
    }

    int totalWedges = secondNeighborsPrefixSum[M];
    // 3: Initialize W to be an array of wedges
    struct Wedge {
        int endpoint1, endpoint2, center;
    };
    vector<Wedge> W(totalWedges);

    //4: parfor u ∈ [0..n):
    for (int u1 = 0; u1 < totalVertices; ++u1) {
        int modDeg = (int)modifiedNeighbors[u1].size();
        // 5: parfor i ← 0 to degu1 (u1) do
        for (int i = 0; i < modDeg; ++i) {
            // We don't need to check rank(v) > rank(u1) because modifiedNeighbors only stores
            // neighbors that are greater rank than the vertex.
            // 6: v ← N (u1)[i] . v = ith neighbor of u1
            int v = modifiedNeighbors[u1][i];
            int base = secondNeighborsPrefixSum[immediateNeighbors[u1] + i];
            int cnt = 0;
            // 7: parfor j ← 0 to degu1 (v) do
            // 8: u2 ← N (v)[j] . u2 = jth neighbor of v (u2 is w)
            for (int w : neighbors[v]) {
                if (ranks[w] > ranks[u1])
                // 9: W [I(i, j)] ← ((u1, u2), 1, v) . (u1, u2) are the endpoints, v is the center of the wedge
                    // u1 and w (u2) are endpoints, v is the center. hence rank(v) > rank(u1)
                    // and rank(w) > rank(u1).
                    W[base + cnt++] = {u1, w, v};
            }
        }
    }

    // ------------- ALGORITHM 3: COUNT-V-WEDGES(W) -------------------------

    /* Parallel sorting (semisort): Collect all wedges into a sequence of key-value pairs with
    key (v, w), then use a parallel sort. Equal keys become contiguous, and each group of size f
    implies (fC2) butterflies.

    NOTE: USE BATCHING FOR PARALLELISM, Partition vertices into fixed-size batches. Each thread takes a batch of
    vertices, enumerates all wedges incident to them, and serially tallies them in a local array or
    map.

    THIS ACHIEVES BETTER PERFORMANCE. SORTING IS O(N log N), THE SLOWEST METHOD FOR WEDGE
    AGGREGATION.
    */
    sort(W.begin(), W.end(), [](auto &a, auto &b){
        if (a.endpoint1 != b.endpoint1) return a.endpoint1 < b.endpoint1;
        return a.endpoint2 < b.endpoint2;
    });

    // Generate endpointWedgeFrequency (count of Wedges on a pair of endpoints)
    // And uniqueWedgeIndexes (these tell us how many wedges there are on UNIQUE endpoints)
    
    vector<pair<pair<int, int>, int>> endpointWedgeFrequency;
    vector<int> uniqueWedgeIndexes;

    /* This is (R, F ) ← GET-FREQ(W ) . Aggregate W and retrieve wedge
    frequencies
    -- Turn this into a function later
    */
    for (int i = 0; i < W.size();) {
        int u1 = W[i].endpoint1;
        int u2 = W[i].endpoint2;
        int j = i + 1;
        
        // Since the endpoints are in contiguous memory space, we scan
        // until they are no longer available
        while (j < W.size() && W[j].endpoint1 == u1 && W[j].
                endpoint2 == u2)
                ++j;
        int freq = j - i;
        // freq is the total count of wedges incident on a specific
        // endpoint pair.

        endpointWedgeFrequency.emplace_back(make_pair(u1, u2), freq);
        uniqueWedgeIndexes.push_back(i);
        i = j;
    }

    vector<long long> butterfliesPerVertex(totalVertices, 0LL);
    long long totalButterflies = 0LL;
    
    for (int i = 0; i < uniqueWedgeIndexes.size(); ++i) {
        /*
            4: parfor i ← 0 to |F | − 1 do
            5: ((u1, u2), d) ← R[i] . u1 and u2 are the wedge endpoints
            6: Store (u1, (dC2)) and (u2, (dC2)) in B. Store butterfly counts per endpoint
        */
        auto freqPair = endpointWedgeFrequency[i];
        int u1 = freqPair.first.first;
        int u2 = freqPair.first.second;
        int freq = freqPair.second;

        // Taking a simple combination instead of adding every single wedge on the same
        // incident endpoints individually.
        long long butterfliesForEndpoints = freq * (freq - 1) / 2;

        butterfliesPerVertex[u1] += butterfliesForEndpoints;
        butterfliesPerVertex[u2] += butterfliesForEndpoints;
        totalButterflies += butterfliesForEndpoints;
        /*
            7: parfor j ← F [i] to F [i + 1] do
            8: ( , , v) ← W [j] . v is the wedge center
            9: Store (v, d − 1) in B
        */
        for (int j = uniqueWedgeIndexes[i]; j < uniqueWedgeIndexes[i + 1]; ++j) {
            int v = W[j].center;
            // Do not need to add this to total butterfly count as that's already included by the
            // vertex counts
            butterfliesPerVertex[v] += freq - 1;
        }
    }

    cout << "Total number of butterflies = "
     << totalButterflies << "\n";
    // End of COUNT-V-WEDGES
    
    return 0;
}