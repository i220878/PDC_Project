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

// Algorithm 5, lines 2–9: UPDATE-V(G = (U, V, E), B, A)
void updateV(
    const vector<vector<int>>& neighbors, // adjacency list for all vertices
    const vector<int>&         ranks,     // vertex ranks
    vector<long long>&         B,         // butterfly counts per vertex (to update)
    const vector<int>&         peelSet,   // A: vertices being peeled this round
    const string&              peelType,  // aggregation method: SORT/HASH/BATCH
    vector<Wedge>&             wedgeBuffer // buffer for temporary wedge storage
) {
    // (line 2) Clear wedge buffer for this round
    wedgeBuffer.clear();

    #pragma omp parallel
    {
        vector<Wedge> localWedges;
        #pragma omp for schedule(dynamic)
        /*
            3: parfor u1 ∈ A do
            4: parfor v ∈ N (u1) do
            5: parfor u2 ∈ N (v) where u2 6 = u1 do
        */
        for (int idx = 0; idx < (int)peelSet.size(); ++idx) {
            int u1 = peelSet[idx];
            for (int v : neighbors[u1]) {
                for (int u2 : neighbors[v]) {
                    if (u2 == u1) continue; // skip self
                    // Only consider wedges where both v and u2 have higher rank than u1
                    if (ranks[v] > ranks[u1] && ranks[u2] > ranks[u1]) {
                        // 6: Store ((u1, u2), 1, v) in W . (u1, u2) is the key, 1 is the frequency
                        localWedges.push_back({u1, u2, v});
                    }
                }
            }
        }
        // Merge local wedges into global buffer
        #pragma omp critical
        wedgeBuffer.insert(wedgeBuffer.end(), localWedges.begin(), localWedges.end());
    }

    // 7: B′ ← COUNT-V-WEDGES(G, W )
    vector<pair<pair<int,int>,int>> wedgeGroups;
    vector<int> groupOffsets;
    if      (peelType == "SORT")  aggregateWedgesSort(wedgeBuffer, wedgeGroups, groupOffsets);
    else if (peelType == "HASH")  aggregateWedgesHash(wedgeBuffer, wedgeGroups, groupOffsets);
    else if (peelType == "BATCH") aggregateWedgesBatch(wedgeBuffer, wedgeGroups, groupOffsets);

    // 8: Subtract corresponding counts B′ from B
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)wedgeGroups.size(); ++i) {
        int u1 = wedgeGroups[i].first.first;
        int u2 = wedgeGroups[i].first.second;
        int freq = wedgeGroups[i].second;
        long long butterflyCount = (long long)freq * (freq - 1) / 2;
        // Subtract from endpoints
        #pragma omp atomic
        B[u1] -= butterflyCount;
        #pragma omp atomic
        B[u2] -= butterflyCount;
        // Subtract from centers
        for (int j = groupOffsets[i]; j < groupOffsets[i+1]; ++j) {
            int v = wedgeBuffer[j].center;
            #pragma omp atomic
            B[v] -= (freq - 1);
        }
    }
}



// Algorithm 5, lines 10–18: PEEL-V(G = (U, V, E), B)
void peelV(
    int                         v1,        // number of vertices in U
    const vector<vector<int>>&  neighbors, // adjacency list
    const vector<int>&          ranks,     // vertex ranks
    vector<long long>&          B,         // butterfly counts per vertex (to update)
    const string&               peelType   // aggregation method: SORT/HASH/BATCH
) {
    // 11: Let K be a bucketing structure mapping U to buckets based on # of butterflies
    long long maxButterflies = 0;
    for (int u = 0; u < v1; ++u)
        maxButterflies = max(maxButterflies, B[u]);
    Bucket K(v1, (int)maxButterflies, true); // extractMin = true

    // Initialize bucket scores with current butterfly counts
    vector<int> initialScores(v1);
    for (int u = 0; u < v1; ++u) initialScores[u] = (int)B[u];
    K.initialize(initialScores);

    // Track which vertices have been peeled
    vector<char> peeled(v1, 0);
    vector<int>  peelSet; // A: vertices to peel in this round
    vector<Wedge> wedgeBuffer;
    // 12: f ← 0
    int numPeeled = 0;

    // 13: while f < | U | do
    while (numPeeled < v1) {
        // 14: A (peelSet) ← all vertices in next bucket (to be peeled)
        peelSet.clear();
        int u;
        while ((u = K.pop()) != -1) {
            if (!peeled[u]) {
                peeled[u] = 1;
                peelSet.push_back(u);
            }
        }
        // 15: f ← f + |A|
        numPeeled += (int)peelSet.size();

        // 16: B ←UPDATE-V(G, B, A) . Update # butterflies
        updateV(neighbors, ranks, B, peelSet, peelType, wedgeBuffer);

        // 17: Update the buckets of changed vertices in B
        #pragma omp parallel for schedule(dynamic)
        for (int idx = 0; idx < (int)peelSet.size(); ++idx) {
            int u1 = peelSet[idx];
            for (int v : neighbors[u1]) {
                for (int u2 : neighbors[v]) {
                    if (u2 < v1 && !peeled[u2]) {
                        long long newScore = B[u2];
                        long long oldScore = K.getScore(u2);
                        int delta = (int)(newScore - oldScore);
                        while (delta > 0) { K.update(u2, +1); --delta; }
                        while (delta < 0) { K.update(u2, -1); ++delta; }
                    }
                }
            }
        }
    }
}

/*
    ALGORITHM 6 EDGE DECOMPOSITION
*/
void updateE(
    const vector<vector<int>>&  neighbors,
    const vector<pair<int,int>>& edgesList,
    const unordered_map<pair<int,int>,int,PairHash>& edgeIdOf,
    const vector<vector<int>>&  incidentEdges,
    vector<long long>&          B_edge,         // per-edge butterfly counts
    const vector<int>&          peelEdges,      // A: list of eids peeled this round
    const string&               peelType,
    vector<pair<pair<int,int>,int>>& contribBuf // temporary buffer of (eid,Δ)
) {
    // (2) Initialize B′ to store updated butterfly counts
    contribBuf.clear();

    // (3) parfor (u1, v1) ∈ A do
    #pragma omp parallel
    {
        vector<pair<pair<int,int>,int>> local;
        #pragma omp for schedule(dynamic)
        for (int idx = 0; idx < (int)peelEdges.size(); ++idx) {
            int eid = peelEdges[idx];
            int u1 = edgesList[eid].first;
            int v1 = edgesList[eid].second;

            // (4) parfor u2 ∈ N(v1) where u2 ≠ u1 do
            for (int u2 : neighbors[v1]) {
                if (u2 == u1) continue;

                // (5) N ← INTERSECT(N(u1), N(u2))
                const auto& Nu1 = neighbors[u1];
                const auto& Nu2 = neighbors[u2];
                vector<int> N;
                int i = 0, j = 0;
                while (i < (int)Nu1.size() && j < (int)Nu2.size()) {
                    if      (Nu1[i] < Nu2[j]) ++i;
                    else if (Nu1[i] > Nu2[j]) ++j;
                    else {
                        N.push_back(Nu1[i]);
                        ++i; ++j;
                    }
                }

                // (6) Store ((u2, v1), |N| - 1) in B′
                int eid2 = edgeIdOf.at({u2, v1});
                local.emplace_back(make_pair(eid2, 1), (int)N.size() - 1);

                // (7) parfor v2 ∈ N where v2 ≠ v1 do
                for (int v2 : N) {
                    if (v2 == v1) continue;
                    // (8) Store ((u1, v2), 1) in B′
                    int eA = edgeIdOf.at({u1, v2});
                    local.emplace_back(make_pair(eA, 1), 1);
                    // (9) Store ((u2, v2), 1) in B′
                    int eB = edgeIdOf.at({u2, v2});
                    local.emplace_back(make_pair(eB, 1), 1);
                }
            }
        }
        // Merge local buffer into global buffer
        #pragma omp critical
        contribBuf.insert(contribBuf.end(), local.begin(), local.end());
    }

    // (10) (B′′, ) ← GET-FREQ(B′)
    vector<pair<pair<int,int>, long long>> R;
    if      (peelType == "SORT") R = aggregateKeyValueSort(contribBuf);
    else if (peelType == "HASH") R = aggregateKeyValueHash(contribBuf);
    else                        R = aggregateKeyValueBatch(contribBuf);

    // (11) Subtract corresponding counts in B′′ from B
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)R.size(); ++i) {
        int eid = R[i].first.first;
        int delta = R[i].second;
        #pragma omp atomic
        B_edge[eid] -= delta;
    }
    // (12) return B (done by reference)
}

/*
    ALGORITHM 6 PEEL-E (Wing Decomposition)
    13: procedure PEEL-E(G = (U, V, E), B)
    14: Let K be a bucketing structure mapping E to buckets based on # of butterflies
    15: f ← 0
    16: while f < m do
    17:   A ← all edges in next bucket (to be peeled)
    18:   f ← f + |A|
    19:   B ← UPDATE-E(G, B, A) . Update # butterflies
    20:   Update the buckets of changed edges in B
    21: return K
*/
void peelE(
    int                                 totalEdges,
    const vector<vector<int>>&          neighbors,
    const vector<pair<int,int>>&        edgesList,
    const unordered_map<pair<int,int>,int,PairHash>& edgeIdOf,
    const vector<vector<int>>&          incidentEdges,
    vector<long long>&                  B_edge,
    const string&                       peelType
) {
    // (14) Build bucket structure for edges based on butterfly counts
    long long maxB = 0;
    for (int i = 0; i < totalEdges; ++i) {
        if (B_edge[i] > maxB) maxB = B_edge[i];
    }
    Bucket K(totalEdges, (int)maxB, /*extractMin=*/true);

    // Initialize bucket scores with current butterfly counts
    vector<int> initScores(totalEdges);
    for (int e = 0; e < totalEdges; ++e) {
        initScores[e] = (int)B_edge[e];
    }
    K.initialize(initScores);

    // Track which edges have been peeled
    vector<char> peeled(totalEdges, 0);
    vector<int>  peelSet; // A: edges to peel in this round
    vector<pair<pair<int,int>,int>> contribBuf;
    // 15: f ← 0
    int peeledCount = 0;

    // 16: while f < m do
    while (peeledCount < totalEdges) {
        peelSet.clear();
        int e;
        // 17: A ← all edges in next bucket (to be peeled)
        while ((e = K.pop()) != -1) {
            if (!peeled[e]) {
                peeled[e] = 1;
                peelSet.push_back(e);
            }
        }
        // 18: f ← f + |A|
        peeledCount += (int)peelSet.size();

        // 19: B ← UPDATE-E(G, B, A) . Update # butterflies
        updateE(neighbors, edgesList, edgeIdOf, incidentEdges,
                B_edge, peelSet, peelType, contribBuf);

        // 20: Update the buckets of changed edges in B
        #pragma omp parallel for schedule(dynamic)
        for (int idx = 0; idx < (int)peelSet.size(); ++idx) {
            int e1 = peelSet[idx];
            int u1 = edgesList[e1].first;
            int v1 = edgesList[e1].second;
            // Any edge incident on u1 or v1 whose score changed needs a bucket update
            for (int s = 0; s < 2; ++s) {
                int side = (s == 0) ? u1 : v1;
                for (int j = 0; j < incidentEdges[side].size(); ++j) {
                    int e2 = incidentEdges[side][j];
                    if (!peeled[e2]) {
                        long long newS = B_edge[e2];
                        long long oldS = K.getScore(e2);
                        int delta = (int)(newS - oldS);
                        while (delta > 0) { K.update(e2, +1); --delta; }
                        while (delta < 0) { K.update(e2, -1); ++delta; }
                    }
                }
            }
        }
    }
    // (21) return K (wing numbers are now in K.score)
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
    int    d           = cmd.getOptionIntValue("-d", 5);
    char*  filename    = cmd.getArgument(0);

    bool ok =
       (countType  == "SORT" || countType == "HASH" || countType == "BATCH")
    && (rankType   == "DEG"  || rankType == "SIDE" || rankType == "ADEG")
    && (peelType   == "NONE" || peelType == "SORT" || peelType == "BATCH" || peelType == "HASH")
    && (perType    == "VERT" || perType == "EDGE")
    && (sparseType == "NONE" || sparseType == "COLOR" || sparseType == "EDGE");
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


    int totalVertices = v1 + v2;
    

    vector<pair<int,int>> edgesList;
    edgesList.reserve(e1 + e2);
    // first, all edges from V‐side input
    for (int u = 0; u < v1; ++u) {
      int start = offsetsV[u], end = (u+1<v1 ? offsetsV[u+1] : e1);
      for (int j = start; j < end; ++j) {
        int w = v1 + edgesV[j];               // global vertex ID on U‐side
        edgesList.emplace_back(u, w);
      }
    }

    // Removing duplicates
    size_t oldE = edgesList.size();
    std::sort(edgesList.begin(), edgesList.end());
    edgesList.erase(
        std::unique(edgesList.begin(), edgesList.end()),
        edgesList.end()
    );
    size_t newE = edgesList.size();
    if (newE < oldE) {
        std::cerr << "Warning: removed "
                << (oldE - newE)
                << " duplicate edges\n";
    }

    double scaleFactor = 1.0;
    if (sparseType == "EDGE") {
        // 1) EDGE sparsification: keep each edge with probability p = 1/d
        double p = 1.0 / d;
        scaleFactor = 1.0 / pow(p, 4);
        srand(time(NULL));
        vector<pair<int, int>> filtered;
        filtered.reserve(edgesList.size());
        for (auto& e : edgesList) {
            double r = (double)rand() / RAND_MAX;
            if (r < p)
                filtered.push_back(e);
        }
        edgesList.swap(filtered);
    }
    else if (sparseType == "COLOR") {
        // 2) COLOR sparsification: assign each vertex a color in [0..d)
        // keep edges whose endpoints match
        // a butterfly survives iff its 3 other vertices match the first: p = 1/d^3
        // so scale = 1/p = d^3
        scaleFactor = pow(d, 3);
    
        srand((unsigned)time(NULL));
        vector<int> color(totalVertices);
        for (int u = 0; u < totalVertices; ++u)
            color[u] = rand() % d;
    
        vector<pair<int,int>> filtered;
        filtered.reserve(edgesList.size());
        for (auto &e : edgesList) {
            if (color[e.first] == color[e.second])
                filtered.push_back(e);
        }
        edgesList.swap(filtered);
    } else {
        // 3) No sparsification, but still scale (scaleFactor = 1)
        scaleFactor = 1.0;
    }

    vector<vector<int>> neighbors(totalVertices);
    t1.start();
    // Rebuild neighbors from possibly sparsified edge list
    for (auto& N : neighbors) N.clear();
    for (auto& e : edgesList) {
        int u = e.first, v = e.second;
        neighbors[u].push_back(v);
        neighbors[v].push_back(u);
    }
    t1.stop();
    t1.reportTotal("Constructing neighbors");

    unordered_map<pair<int,int>,int,PairHash> edgeIdOf;
    edgeIdOf.reserve(edgesList.size());
    vector<vector<int>> incidentEdges(totalVertices);  
    for (int eid = 0; eid < (int)edgesList.size(); ++eid) {
      auto [u,v] = edgesList[eid];
      edgeIdOf[{u,v}] = eid;
      // if undirected, also map (v,u) → eid
      edgeIdOf[{v,u}] = eid;
      // record that this eid touches both endpoints
      incidentEdges[u].push_back(eid);
      incidentEdges[v].push_back(eid);
    }

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
        cout << "Total number of butterflies = " << (long long)(totalButterflies * scaleFactor) << "\n";
        
        if (peelType != "NONE") {
            t1.start();
            peelV(v1, neighbors, ranks, B, peelType);
            t1.stop();
            t1.reportTotal("Time To Peel");
        }
        
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
        //*/
        cout << "Total number of butterflies = " << (long long)((totalButterflies / 4) * scaleFactor) << "\n";

        if (peelType != "NONE") {
            t1.start();
        
            // Build per-edge butterfly count array, indexed by edge ID
            // edgeCounts is a vector of ((u,v), count) pairs
            vector<long long> B_edge(edgesList.size(), 0);
            for (const auto& ec : edgeCounts) {
                int eid;
                // Try both (u,v) and (v,u) in case of undirected edges
                auto it = edgeIdOf.find(ec.first);
                if (it != edgeIdOf.end()) {
                    eid = it->second;
                }
                else {
                    auto it2 = edgeIdOf.find({ec.first.second, ec.first.first});
                    if (it2 != edgeIdOf.end())
                        eid = it2->second;
                    else continue;// If not found, skip this edge
                }
                B_edge[eid] = ec.second;
            }
        
            // Call edge peeling routine (Algorithm 6)
            // This will update B_edge in place with wing numbers
            peelE(
                (int)edgesList.size(),   // total number of edges
                neighbors,               // adjacency list
                edgesList,               // edge list (eid -> (u,v))
                edgeIdOf,                // map from (u,v) to eid
                incidentEdges,           // for each vertex, list of incident eids
                B_edge,                  // per-edge butterfly counts (to update)
                peelType                 // aggregation method
            );
        
            t1.stop();
            t1.reportTotal("Time To Peel Edges");
        }
    
        startTime.stop();
        startTime.reportTotal("Total time to run");
    }
    

    return 0;
}