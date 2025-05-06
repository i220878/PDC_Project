#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <cmath>
#include "gettime.h"
#include <omp.h>
#include "parseCommandLine.h"
#include "bucket.h"

using namespace std;

template<typename T>
void printOffsetsAndEdges(std::vector<T> off1,
                          std::vector<T> edge1,
                          std::vector<T> off2,
                          std::vector<T> edge2) {
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


/* RANKING VIA DEGREE ORDER */
void computeDegreeOrder(
    vector<vector<int>>& neighbors,
    vector<int>&            degrees,
    vector<int>&            sortedVertices,
    vector<int>&            ranks,
    vector<vector<int>>&    modifiedNeighbors
) {
    int n = (int)neighbors.size();
    for (int u = 0; u < n; ++u) {
        degrees[u] = (int)neighbors[u].size();
    }
    for (int u = 0; u < n; ++u) {
        sortedVertices[u] = u;
    }
    // Vertices ordered in descending order of degree
    sort(sortedVertices.begin(), sortedVertices.end(),
         [&degrees](int a, int b) {
             return degrees[a] > degrees[b];
         });
    // Ranks is the inverse of sorted vertices
    for (int i = 0; i < n; ++i) {
        ranks[sortedVertices[i]] = i;
    }
    // Neighbors are sorted in descending order of rank
    #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < n; ++u) {
        auto& N = neighbors[u];
        sort(N.begin(), N.end(),
             [&ranks](int a, int b) {
                 return ranks[a] > ranks[b];
             });
    }
    // Modified neighbors has neighbors with rank[neighbor] > vertex,
    // all of a vertex's neighbors will have a rank greater than it 
    #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < n; ++u) {
        auto& M = modifiedNeighbors[u];
        M.clear();
        for (int v : neighbors[u]) {
            if (ranks[v] > ranks[u]) {
                M.push_back(v);
            }
        }
    }
}

/* RANKING VIA APPROX DEGREE ORDER */
void computeApproxDegreeOrder(
    vector<vector<int>>& neighbors,
    vector<int>&            degrees,
    vector<int>&            sortedVertices,
    vector<int>&            ranks,
    vector<vector<int>>&    modifiedNeighbors
) {
    int n = (int)neighbors.size();
    for (int u = 0; u < n; ++u) {
        degrees[u] = (int)neighbors[u].size();
    }
    for (int u = 0; u < n; ++u) {
        sortedVertices[u] = u;
    }
    // Vertices ordered in descending order of degree
    sort(sortedVertices.begin(), sortedVertices.end(),
    [&degrees](int a, int b) {
        int da = floor(log2(degrees[a] > 0 ? degrees[a] : 1));
        int db = floor(log2(degrees[b] > 0 ? degrees[b] : 1));
        if (da != db)
            return da > db;
        return a < b;  // tie‐break on index
    });

    for (int i = 0; i < n; ++i) {
        ranks[sortedVertices[i]] = i;
    }
    // Neighbors are sorted in descending order of rank
    #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < n; ++u) {
        auto& N = neighbors[u];
        sort(N.begin(), N.end(),
             [&ranks](int a, int b) {
                 return ranks[a] > ranks[b];
             });
    }
    // Modified neighbors has neighbors with rank[neighbor] > vertex,
    // all of a vertex's neighbors will have a rank greater than it 
    #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < n; ++u) {
        auto& M = modifiedNeighbors[u];
        M.clear();
        for (int v : neighbors[u]) {
            if (ranks[v] > ranks[u]) {
                M.push_back(v);
            }
        }
    }
}

/* RANKING VIA SIDE ORDER */
void computeSideOrder(
    int v1, int v2,
    const vector<int>& offsetsV, int e1,
    const vector<int>& offsetsU, int e2,
    vector<vector<int>>& neighbors,
    vector<int>&            ranks,
    vector<vector<int>>&    modifiedNeighbors
) {
    int n = v1 + v2;

    // Count wedges in V
    double wedgeV = 0;
    for (int i = 0; i < v1; ++i) {
        int start = offsetsV[i];
        int end = (i + 1 < v1 ? offsetsV[i + 1] : e1);
        int d = end - start;
        wedgeV += d * (d - 1) / 2;
    }

    // Count wedges in U
    double wedgeU = 0;
    for (int j = 0; j < v2; ++j) {
        int start = offsetsU[j];
        int end = (j + 1 < v2 ? offsetsU[j + 1] : e2);
        int d = end - start;
        wedgeU += d * (d - 1) / 2;
    }

    // Assign global ranks
    ranks.resize(n);
    if (wedgeV <= wedgeU) {
        // V side first: V[0..v1-1], then U[v1..v1+v2-1]
        for (int i = 0; i < v1; ++i)
            ranks[i] = i;
        for (int j = 0; j < v2; ++j)
            ranks[v1 + j] = v1 + j;
    } else {
        // U side first: U[v1..v1+v2-1] get ranks [0..v2-1], then V[0..v1-1] get [v2..v2+v1-1]
        for (int j = 0; j < v2; ++j)
            ranks[v1 + j] = j;
        for (int i = 0; i < v1; ++i)
            ranks[i] = v2 + i;
    }

    // Neighbors are sorted in descending order 1of rank
    #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < n; ++u) {
        auto& N = neighbors[u];
        sort(N.begin(), N.end(),
             [&ranks](int a, int b) {
                 return ranks[a] > ranks[b];
             });
    }
    // Modified neighbors has neighbors with rank[neighbor] > vertex,
    // all of a vertex's neighbors will have a rank greater than it 
    modifiedNeighbors.resize(n);
    #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < n; ++u) {
        auto& M = modifiedNeighbors[u];
        M.clear();
        for (int v : neighbors[u]) {
            if (ranks[v] > ranks[u]) {
                M.push_back(v);
            }
        }
    }
}

struct Wedge {
    int endpoint1, endpoint2, center;
};

struct PairHash {
    size_t operator()(pair<int,int> const &p) const {
      return (size_t)p.first * 1000003u + p.second;
    }
};

// Aggregate wedges by (u1,u2) using sorting
void aggregateWedgesSort(
    const vector<Wedge>& wedges,
    vector<pair<pair<int,int>,int>>& groups,
    vector<int>& groupOffsets
) {
    int m = (int)wedges.size();
    vector<Wedge> sortedWedges = wedges;
    sort(sortedWedges.begin(), sortedWedges.end(),
         [](const Wedge& a, const Wedge& b) {
             if (a.endpoint1 != b.endpoint1) return a.endpoint1 < b.endpoint1;
             return a.endpoint2 < b.endpoint2;
         });
    groups.clear();
    groupOffsets.clear();
    groupOffsets.push_back(0);
    for (int i = 0; i < m; ) {
        int u1 = sortedWedges[i].endpoint1, u2 = sortedWedges[i].endpoint2;
        int j = i + 1;
        while (j < m && sortedWedges[j].endpoint1 == u1 && sortedWedges[j].endpoint2 == u2) ++j;
        int freq = j - i;
        groups.emplace_back(make_pair(u1, u2), freq);
        groupOffsets.push_back(j);
        i = j;
    }
}

// Aggregate wedges by (u1,u2) using parallel hash maps
void aggregateWedgesHash(
    const vector<Wedge>& wedges,
    vector<pair<pair<int,int>,int>>& groups,
    vector<int>& groupOffsets
) {
    int m = (int)wedges.size();
    int T = omp_get_max_threads();
    vector<unordered_map<pair<int,int>, int, PairHash>> localMaps(T);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& myMap = localMaps[tid];
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < m; ++i) {
            auto key = make_pair(wedges[i].endpoint1, wedges[i].endpoint2);
            myMap[key]++;
        }
    }
    unordered_map<pair<int,int>, int, PairHash> freq;
    freq.reserve(m);
    for (int t = 0; t < T; ++t) {
        for (auto& kv : localMaps[t]) {
            freq[kv.first] += kv.second;
        }
    }
    groups.clear();
    groupOffsets.clear();
    groupOffsets.push_back(0);
    vector<pair<int,int>> keys;
    keys.reserve(freq.size());
    for (auto& kv : freq) keys.push_back(kv.first);
    sort(keys.begin(), keys.end());
    int pos = 0;
    for (auto& k : keys) {
        int d = freq[k];
        groups.emplace_back(k, d);
        pos += d;
        groupOffsets.push_back(pos);
    }
}

// Aggregate wedges by (u1,u2) using batching and hash maps
void aggregateWedgesBatch(
    const vector<Wedge>& wedges,
    vector<pair<pair<int,int>,int>>& groups,
    vector<int>& groupOffsets
) {
    int m = (int)wedges.size();
    const int B = 1024;
    int T = omp_get_max_threads();
    vector<unordered_map<pair<int,int>, int, PairHash>> localMaps(T);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& myMap = localMaps[tid];
        #pragma omp for schedule(dynamic)
        for (int s = 0; s < m; s += B) {
            int e = min(m, s + B);
            for (int i = s; i < e; ++i) {
                auto key = make_pair(wedges[i].endpoint1, wedges[i].endpoint2);
                myMap[key]++;
            }
        }
    }
    unordered_map<pair<int,int>, int, PairHash> freq;
    freq.reserve(m);
    for (int t = 0; t < T; ++t) {
        for (auto& kv : localMaps[t]) {
            freq[kv.first] += kv.second;
        }
    }
    groups.clear();
    groupOffsets.clear();
    groupOffsets.push_back(0);
    vector<pair<int,int>> keys;
    keys.reserve(freq.size());
    for (auto& kv : freq) keys.push_back(kv.first);
    sort(keys.begin(), keys.end());
    int pos = 0;
    for (auto& k : keys) {
        int d = freq[k];
        groups.emplace_back(k, d);
        pos += d;
        groupOffsets.push_back(pos);
    }
}

// Aggregate key-value pairs by key using sorting
vector<pair<pair<int,int>, long long>>
aggregateKeyValueSort(
    vector<pair<pair<int,int>,int>>& pairs
) {
    sort(pairs.begin(), pairs.end(),
         [](const pair<pair<int,int>,int>& a, const pair<pair<int,int>,int>& b) {
             if (a.first.first != b.first.first) return a.first.first < b.first.first;
             return a.first.second < b.first.second;
         });
    vector<pair<pair<int,int>, long long>> result;
    result.reserve(pairs.size());
    for (int i = 0, n = (int)pairs.size(); i < n; ) {
        auto key = pairs[i].first;
        long long sum = 0;
        int j = i;
        while (j < n && pairs[j].first == key) {
            sum += pairs[j].second;
            ++j;
        }
        result.emplace_back(key, sum);
        i = j;
    }
    return result;
}

// Aggregate key-value pairs by key using parallel hash maps
vector<pair<pair<int,int>, long long>>
aggregateKeyValueHash(
    const vector<pair<pair<int,int>,int>>& pairs
) {
    int m = (int)pairs.size();
    int T = omp_get_max_threads();
    vector<unordered_map<pair<int,int>, long long, PairHash>> localMaps(T);
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < m; ++i) {
        int tid = omp_get_thread_num();
        localMaps[tid][pairs[i].first] += pairs[i].second;
    }
    unordered_map<pair<int,int>, long long, PairHash> freq;
    freq.reserve(m);
    for (int t = 0; t < T; ++t) {
        for (auto& kv : localMaps[t]) {
            freq[kv.first] += kv.second;
        }
    }
    vector<pair<pair<int,int>, long long>> result;
    result.reserve(freq.size());
    vector<pair<int,int>> keys;
    keys.reserve(freq.size());
    for (auto& kv : freq) keys.push_back(kv.first);
    sort(keys.begin(), keys.end());
    for (auto& k : keys) {
        result.emplace_back(k, freq[k]);
    }
    return result;
}

// Aggregate key-value pairs by key using batching and hash maps
vector<pair<pair<int,int>, long long>>
aggregateKeyValueBatch(
    const vector<pair<pair<int,int>,int>>& pairs
) {
    int m = (int)pairs.size();
    const int B = 1024;
    int T = omp_get_max_threads();
    vector<unordered_map<pair<int,int>, long long, PairHash>> localMaps(T);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& myMap = localMaps[tid];
        #pragma omp for schedule(dynamic)
        for (int s = 0; s < m; s += B) {
            int e = min(m, s + B);
            for (int i = s; i < e; ++i) {
                myMap[pairs[i].first] += pairs[i].second;
            }
        }
    }
    unordered_map<pair<int,int>, long long, PairHash> freq;
    freq.reserve(m);
    for (int t = 0; t < T; ++t) {
        for (auto& kv : localMaps[t]) {
            freq[kv.first] += kv.second;
        }
    }
    vector<pair<pair<int,int>, long long>> result;
    result.reserve(freq.size());
    vector<pair<int,int>> keys;
    keys.reserve(freq.size());
    for (auto& kv : freq) keys.push_back(kv.first);
    sort(keys.begin(), keys.end());
    for (auto& k : keys) {
        result.emplace_back(k, freq[k]);
    }
    return result;
}


int main(int argc, char** argv) {
    // parse command‐line
    commandLine cmd(argc, argv,
      "-countType SORT "
      "-rankType DEG "
      "-peelType NONE "
      "-per VERT "
      "-sparseType NONE "
      "-d <int> "
      "<input-file>"
    );
    string countType   = cmd.getOptionValue("-countType",   "BATCH");
    string rankType    = cmd.getOptionValue("-rankType",    "ADEG");
    string peelType    = cmd.getOptionValue("-peelType",    "NONE");
    string perType     = cmd.getOptionValue("-per",         "VERT");
    string sparseType  = cmd.getOptionValue("-sparseType",  "NONE");
    int    d           = cmd.getOptionIntValue("-d", 25);
    char*  filename    = cmd.getArgument(0);

    bool ok =
       (countType  == "SORT" || countType == "HASH" || countType == "BATCH")
    && (rankType   == "DEG"  || rankType == "SIDE" || rankType == "ADEG")
    && (peelType   == "NONE")
    && (perType    == "VERT" || perType == "EDGE")
    && (sparseType == "NONE");
    if (!ok) cmd.badArgument();

    timer startTime;
    startTime.start();
    timer t1;
    t1.start();

    std::ifstream file(filename);
    // Line is first index of file
    std::string line;

    std::vector<int> offsetsV, edgesV;
    std::vector<int> offsetsU, edgesU;

    // Iterate over first 5 linse of the input, containing "AdjacencyHypergraph",
    // nv, mv, nu, mu
    // v1 = Number of Vertices for bipartition V,
    // e1 = Number of Edges for bipartition V,
    // v2 = Number of Vertices for bipartition U,
    // e2 = Number of Edges for bipartition U
    int v1,e1,v2,e2;
    getline(file, line);
    if (line != "AdjacencyHypergraph") {
        std::cerr << "Error: First line must be 'AdjacencyHypergraph'"
                  << std::endl;
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

    t1.stop();
    t1.reportTotal("file_reading");

    t1.start();

    int totalVertices = v1 + v2;
    vector<vector<int>> neighbors(totalVertices);
    // All neighbros are given by the edges from a node in V to a node in U
    for (int i = 0; i < v1; ++i) {
        int start = offsetsV[i];
        int end   = (i + 1 < v1 ? offsetsV[i + 1] : e1);
        for (int j = start; j < end; ++j) {
            int edgeId = edgesV[j];      // V's neighbor in U
            neighbors[i].push_back(v1 + edgeId);
            neighbors[v1 + edgeId].push_back(i); // Since the graph is undirected, edges go both ways
        }
    }
    t1.stop();
    t1.reportTotal("Constructing neighbors");

    /*
        STEP 1: COMPUTE VERTEX RANK
    */
    t1.start();
    vector<int> degrees(totalVertices), sortedVertices(totalVertices), ranks(totalVertices);
    vector<vector<int>> modifiedNeighbors(totalVertices);

    if (rankType == "DEG") {
        computeDegreeOrder(neighbors, degrees, sortedVertices, ranks, modifiedNeighbors);
        }
    else if (rankType == "SIDE") {
        computeSideOrder(v1, v2, offsetsV, e1, offsetsU, e2, neighbors, ranks, modifiedNeighbors);
    }
    else if (rankType == "ADEG") {
        computeApproxDegreeOrder(neighbors, degrees, sortedVertices, ranks, modifiedNeighbors);
    }
    else {
        cmd.badArgument();
    }
    t1.stop();
    t1.reportTotal("Sorting and modified neighborhood shenanigans");
    /*
        We define a modified degree,
        degv (u), to be the number of neighbors u′ ∈ N (u) such that
        rank(u′) > rank(v)

        NOTE: MODIFIED DEGREE is simply modifiedNeighbors[index].size()
    */

    /*
        STEP 2: WEDGE RETRIEVAL / ENUMERATION
    */
    // Build modifiedNeighbors based on ranks
    // 2: Use PREFIX-SUM to compute a function I that maps wedges to indices in order
    // Build immediateNeighbors: exclusive prefix‐sum of modified‐degrees
    t1.start();
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
    #pragma omp parallel for schedule(dynamic)
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
    t1.stop();
    t1.reportTotal("Calculating prefix sums");

    t1.start();
    int totalWedges = (int)secondNeighborsPrefixSum[M];
    // 3: Initialize W to be an array of wedges
    vector<Wedge> W(totalWedges);
    t1.stop();
    t1.reportTotal("Initialising Wedge array");

    t1.start();
    //4: parfor u ∈ [0..n):
    
    /*\ Debugging I did to check if different rankTypes were effecting the creation
    of modifiedNeighbors
    cout << "\nOUTPUTTING MODIFIED NEIGHBORS SIZE: " << modifiedNeighbors.size();
    int sum = 0;
    for (int i = 0; i < modifiedNeighbors.size(); ++i) {
        sum += modifiedNeighbors[i].size();
    }
    cout << "\nTOTAL NEIGHBORS IN EACH MODIFIED NEIGHBOR: " << sum << '\n';
    //*/
    #pragma omp parallel for schedule(dynamic)
    for (int u1 = 0; u1 < totalVertices; ++u1) {
        int modDeg = modifiedNeighbors[u1].size();
        // 5: parfor i ← 0 to degu1 (u1) do
        for (int i = 0; i < modDeg; ++i) {
            // We don't need to check rank(v) > rank(u1) because modifiedNeighbors only stores
            // neighbors that are greater rank than the vertex.
            // 6: v ← N (u1)[i] . v = ith neighbor of u1
            int v = modifiedNeighbors[u1][i];
            int base = (int)secondNeighborsPrefixSum[
                          immediateNeighbors[u1] + i];
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
    t1.stop();
    t1.reportTotal("Wedge retrieval");

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
    // Algorithm 3: Parallel work‐efficient butterfly counting per‐vertex
    if (perType == "VERT") {
        // 2: (R, F) ← GET‐FREQ(W) . Aggregate W and retrieve wedge frequencies
        vector<pair<pair<int,int>,int>> R;
        vector<int>                     F;
        if      (countType == "SORT")  aggregateWedgesSort(W, R, F);
        else if (countType == "HASH")  aggregateWedgesHash(W, R, F);
        else if (countType == "BATCH") aggregateWedgesBatch(W, R, F);
        else cmd.badArgument();
    
        // 3: Initialize B to store butterfly counts per vertex
        vector<long long> B(totalVertices, 0LL);
        long long         totalButterflies = 0;
    
        // 4: parfor i = 0 to |R|−1 do
        // 5:   ((u1, u2), d) ← R[i] . u1 and u2 are the wedge endpoint
        // 6:   Store (u1, (d choose 2)) and (u2, (d choose 2)) in B . Store butterfly counts per endpoint
        #pragma omp parallel for reduction(+ : totalButterflies)
        for (int i = 0; i < (int)R.size(); ++i) {
            int u1 = R[i].first.first;
            int u2 = R[i].first.second;
            int d  = R[i].second;
            long long c = (long long)d * (d - 1) / 2;
            #pragma omp atomic
            B[u1] += c;
            #pragma omp atomic
            B[u2] += c;
            totalButterflies += c;
        }
    
    
        #pragma omp parallel for
        for (int i = 0; i < (int)R.size(); ++i) {
            int start = F[i], end = F[i+1];
            int d     = R[i].second;
            long long c = d - 1;
            // 7: parfor j ← F [i] to F [i + 1] do
            // 8:   ( , , v) ← W [j] . v is the wedge center
            for (int j = start; j < end; ++j) {
                int v = W[j].center;
                #pragma omp atomic
                //9: Store (v, d − 1) in B . Store butterfly counts per center
                B[v] += c;
            }
        }
        // (B, ) ← GET-FREQ(B) . Aggregate B and get butterfly counts
    
        t1.stop();
        t1.reportTotal("Calculating per‐vertex butterfly counts");
        cout << "Total number of butterflies = " << totalButterflies << "\n";
        startTime.stop();
        startTime.reportTotal("Total time to run");
    }
    // Algorithm 4: Parallel work‐efficient butterfly counting per‐edge
    else if (perType == "EDGE") {
        // 2: (R, F) ← GET‐FREQ(W) . Aggregate wedges W into groups by (u1,u2)
        vector<pair<pair<int,int>,int>> R;
        vector<int>                     F;
        if      (countType == "SORT")  aggregateWedgesSort(W, R, F);
        else if (countType == "HASH")  aggregateWedgesHash(W, R, F);
        else if (countType == "BATCH") aggregateWedgesBatch(W, R, F);
        else cmd.badArgument();
    
        int T = omp_get_max_threads();
        vector<vector<pair<pair<int,int>,int>>> local_B_contribs(T);
    
        size_t totalContribs = 0;
        for (int i = 0; i < (int)R.size(); ++i) {
            totalContribs += 2 * (F[i+1] - F[i]); // 2 contributions per wedge
        }
        
        // Pre-allocate for each thread
        for (int t = 0; t < T; ++t) {
            local_B_contribs[t].reserve(totalContribs / T + 1);
        }
    
        // 4: parfor i = 0 to |R|−1 do
        // 5:  ((u1, u2), d) ← R[i]
        // 6:   for j = F[i] to F[i+1]−1:
        // 7:     (u1,u2, v) ← W[j], store ((u1,v),(d−1)) and ((u2,v),(d−1))
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto& my_contribs = local_B_contribs[tid];
            
            #pragma omp for schedule(dynamic)
            for (int i = 0; i < (int)R.size(); ++i) {
                int u1    = R[i].first.first;
                int u2    = R[i].first.second;
                int d     = R[i].second;
                int start = F[i], end = F[i+1];
                int val   = d - 1;
                for (int j = start; j < end; ++j) {
                    int v = W[j].center;
                    my_contribs.emplace_back(make_pair(u1, v), val);
                    my_contribs.emplace_back(make_pair(u2, v), val);
                }
            }
        }
    
        // Merge thread-local vectors (serially to avoid races)
        vector<pair<pair<int,int>,int>> B_contribs;
        B_contribs.reserve(totalContribs);
        for (int t = 0; t < T; ++t) {
            B_contribs.insert(B_contribs.end(), 
                            local_B_contribs[t].begin(), 
                            local_B_contribs[t].end());
        }

        // 9: (B_edge, _) ← GET‐FREQ(B_contribs)
        vector<pair<pair<int,int>, long long>> edgeCounts;
        if      (countType == "SORT")  edgeCounts = aggregateKeyValueSort(B_contribs);
        else if (countType == "HASH")  edgeCounts = aggregateKeyValueHash(B_contribs);
        else if (countType == "BATCH") edgeCounts = aggregateKeyValueBatch(B_contribs);

        // return B
        long long totalButterflies = 0;
        for (auto &ec : edgeCounts) {
            totalButterflies += ec.second;
        }
        cout << "Total number of butterflies = " << totalButterflies / 4 << "\n";
    
        startTime.stop();
        startTime.reportTotal("Total time to run");
    }
    

    return 0;
}
