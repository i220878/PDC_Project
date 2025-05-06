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
// #include <omp.h>
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
    int w = 2; // adjust if your n > 99

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
    cout << string(w+1, ' ');
    for (int j = 0; j < n; ++j) {
        cout << string(w, '-') << '-';
        if (j == v1-1) cout << "-";
    }
    cout << "\n";

    recount = false;
    for (int i = 0; i < n; ++i) {
        if (recount) cout << setw(w) << i - v1 << ' ';
        else cout << setw(w) << i << ' ';
        cout << "| ";
        if (i == v1 - 1) {
            recount = true;
        }
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
    sort(sortedVertices.begin(), sortedVertices.end(),
         [&degrees](int a, int b) {
             return degrees[a] > degrees[b];
         });
    for (int i = 0; i < n; ++i) {
        ranks[sortedVertices[i]] = i;
    }
    // #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < n; ++u) {
        auto& N = neighbors[u];
        sort(N.begin(), N.end(),
             [&ranks](int a, int b) {
                 return ranks[a] > ranks[b];
             });
    }
    // #pragma omp parallel for schedule(dynamic)
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
    sort(sortedVertices.begin(), sortedVertices.end(),
    [&degrees](int a, int b) {
        int da = floor(log2(degrees[a] > 0 ? degrees[a] : 1));
        int db = floor(log2(degrees[b] > 0 ? degrees[b] : 1));
        if (da != db)
            return da > db;
        return a < b;
    });

    for (int i = 0; i < n; ++i) {
        ranks[sortedVertices[i]] = i;
    }
    // #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < n; ++u) {
        auto& N = neighbors[u];
        sort(N.begin(), N.end(),
             [&ranks](int a, int b) {
                 return ranks[a] > ranks[b];
             });
    }
    // #pragma omp parallel for schedule(dynamic)
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

    double wedgeV = 0;
    for (int i = 0; i < v1; ++i) {
        int start = offsetsV[i];
        int end = (i + 1 < v1 ? offsetsV[i + 1] : e1);
        int d = end - start;
        wedgeV += d * (d - 1) / 2;
    }

    double wedgeU = 0;
    for (int j = 0; j < v2; ++j) {
        int start = offsetsU[j];
        int end = (j + 1 < v2 ? offsetsU[j + 1] : e2);
        int d = end - start;
        wedgeU += d * (d - 1) / 2;
    }

    ranks.resize(n);
    if (wedgeV <= wedgeU) {
        for (int i = 0; i < v1; ++i)
            ranks[i] = i;
        for (int j = 0; j < v2; ++j)
            ranks[v1 + j] = v1 + j;
    } else {
        for (int j = 0; j < v2; ++j)
            ranks[v1 + j] = j;
        for (int i = 0; i < v1; ++i)
            ranks[i] = v2 + i;
    }

    // #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < n; ++u) {
        auto& N = neighbors[u];
        sort(N.begin(), N.end(),
             [&ranks](int a, int b) {
                 return ranks[a] > ranks[b];
             });
    }
    modifiedNeighbors.resize(n);
    // #pragma omp parallel for schedule(dynamic)
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

void aggregateWedgesHash(
    const vector<Wedge>& wedges,
    vector<pair<pair<int,int>,int>>& groups,
    vector<int>& groupOffsets
) {
    int m = (int)wedges.size();
    int T = 1; // omp_get_max_threads();
    vector<unordered_map<pair<int,int>, int, PairHash>> localMaps(T);
    // #pragma omp parallel
    {
        int tid = 0; // omp_get_thread_num();
        auto& myMap = localMaps[tid];
        // #pragma omp for schedule(dynamic)
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

void aggregateWedgesBatch(
    const vector<Wedge>& wedges,
    vector<pair<pair<int,int>,int>>& groups,
    vector<int>& groupOffsets
) {
    int m = (int)wedges.size();
    const int B = 1024;
    int T = 1; // omp_get_max_threads();
    vector<unordered_map<pair<int,int>, int, PairHash>> localMaps(T);
    // #pragma omp parallel
    {
        int tid = 0; // omp_get_thread_num();
        auto& myMap = localMaps[tid];
        // #pragma omp for schedule(dynamic)
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

vector<pair<pair<int,int>, long long>>
aggregateKeyValueHash(
    const vector<pair<pair<int,int>,int>>& pairs
) {
    int m = (int)pairs.size();
    int T = 1; // omp_get_max_threads();
    vector<unordered_map<pair<int,int>, long long, PairHash>> localMaps(T);
    // #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < m; ++i) {
        int tid = 0; // omp_get_thread_num();
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

vector<pair<pair<int,int>, long long>>
aggregateKeyValueBatch(
    const vector<pair<pair<int,int>,int>>& pairs
) {
    int m = (int)pairs.size();
    const int B = 1024;
    int T = 1; // omp_get_max_threads();
    vector<unordered_map<pair<int,int>, long long, PairHash>> localMaps(T);
    // #pragma omp parallel
    {
        int tid = 0; // omp_get_thread_num();
        auto& myMap = localMaps[tid];
        // #pragma omp for schedule(dynamic)
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

void updateV(
    const vector<vector<int>>& neighbors,
    const vector<int>&         ranks,
    vector<long long>&         B,
    const vector<int>&         peelSet,
    const string&              peelType,
    vector<Wedge>&             wedgeBuffer
) {
    wedgeBuffer.clear();

    // #pragma omp parallel
    {
        vector<Wedge> localWedges;
        // #pragma omp for schedule(dynamic)
        for (int idx = 0; idx < (int)peelSet.size(); ++idx) {
            int u1 = peelSet[idx];
            for (int v : neighbors[u1]) {
                for (int u2 : neighbors[v]) {
                    if (u2 == u1) continue;
                    if (ranks[v] > ranks[u1] && ranks[u2] > ranks[u1]) {
                        localWedges.push_back({u1, u2, v});
                    }
                }
            }
        }
        // #pragma omp critical
        wedgeBuffer.insert(wedgeBuffer.end(), localWedges.begin(), localWedges.end());
    }

    vector<pair<pair<int,int>,int>> wedgeGroups;
    vector<int> groupOffsets;
    if      (peelType == "SORT")  aggregateWedgesSort(wedgeBuffer, wedgeGroups, groupOffsets);
    else if (peelType == "HASH")  aggregateWedgesHash(wedgeBuffer, wedgeGroups, groupOffsets);
    else if (peelType == "BATCH") aggregateWedgesBatch(wedgeBuffer, wedgeGroups, groupOffsets);

    // #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)wedgeGroups.size(); ++i) {
        int u1 = wedgeGroups[i].first.first;
        int u2 = wedgeGroups[i].first.second;
        int freq = wedgeGroups[i].second;
        long long butterflyCount = (long long)freq * (freq - 1) / 2;
        // #pragma omp atomic
        B[u1] -= butterflyCount;
        // #pragma omp atomic
        B[u2] -= butterflyCount;
        for (int j = groupOffsets[i]; j < groupOffsets[i+1]; ++j) {
            int v = wedgeBuffer[j].center;
            // #pragma omp atomic
            B[v] -= (freq - 1);
        }
    }
}

void peelV(
    int                         v1,
    const vector<vector<int>>&  neighbors,
    const vector<int>&          ranks,
    vector<long long>&          B,
    const string&               peelType
) {
    long long maxButterflies = 0;
    for (int u = 0; u < v1; ++u)
        maxButterflies = max(maxButterflies, B[u]);
    Bucket K(v1, (int)maxButterflies, true);

    vector<int> initialScores(v1);
    for (int u = 0; u < v1; ++u) initialScores[u] = (int)B[u];
    K.initialize(initialScores);

    vector<char> peeled(v1, 0);
    vector<int>  peelSet;
    vector<Wedge> wedgeBuffer;
    int numPeeled = 0;

    while (numPeeled < v1) {
        peelSet.clear();
        int u;
        while ((u = K.pop()) != -1) {
            if (!peeled[u]) {
                peeled[u] = 1;
                peelSet.push_back(u);
            }
        }
        numPeeled += (int)peelSet.size();

        updateV(neighbors, ranks, B, peelSet, peelType, wedgeBuffer);

        // #pragma omp parallel for schedule(dynamic)
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

void updateE(
    const vector<vector<int>>&  neighbors,
    const vector<pair<int,int>>& edgesList,
    const unordered_map<pair<int,int>,int,PairHash>& edgeIdOf,
    const vector<vector<int>>&  incidentEdges,
    vector<long long>&          B_edge,
    const vector<int>&          peelEdges,
    const string&               peelType,
    vector<pair<pair<int,int>,int>>& contribBuf
) {
    contribBuf.clear();

    // #pragma omp parallel
    {
        vector<pair<pair<int,int>,int>> local;
        // #pragma omp for schedule(dynamic)
        for (int idx = 0; idx < (int)peelEdges.size(); ++idx) {
            int eid = peelEdges[idx];
            int u1 = edgesList[eid].first;
            int v1 = edgesList[eid].second;

            for (int u2 : neighbors[v1]) {
                if (u2 == u1) continue;

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

                int eid2 = edgeIdOf.at({u2, v1});
                local.emplace_back(make_pair(eid2, 1), (int)N.size() - 1);

                for (int v2 : N) {
                    if (v2 == v1) continue;
                    int eA = edgeIdOf.at({u1, v2});
                    local.emplace_back(make_pair(eA, 1), 1);
                    int eB = edgeIdOf.at({u2, v2});
                    local.emplace_back(make_pair(eB, 1), 1);
                }
            }
        }
        // #pragma omp critical
        contribBuf.insert(contribBuf.end(), local.begin(), local.end());
    }

    vector<pair<pair<int,int>, long long>> R;
    if      (peelType == "SORT") R = aggregateKeyValueSort(contribBuf);
    else if (peelType == "HASH") R = aggregateKeyValueHash(contribBuf);
    else                        R = aggregateKeyValueBatch(contribBuf);

    // #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)R.size(); ++i) {
        int eid = R[i].first.first;
        int delta = R[i].second;
        // #pragma omp atomic
        B_edge[eid] -= delta;
    }
}

void peelE(
    int                                 totalEdges,
    const vector<vector<int>>&          neighbors,
    const vector<pair<int,int>>&        edgesList,
    const unordered_map<pair<int,int>,int,PairHash>& edgeIdOf,
    const vector<vector<int>>&          incidentEdges,
    vector<long long>&                  B_edge,
    const string&                       peelType
) {
    long long maxB = 0;
    for (int i = 0; i < totalEdges; ++i) {
        if (B_edge[i] > maxB) maxB = B_edge[i];
    }
    Bucket K(totalEdges, (int)maxB, /*extractMin=*/true);

    vector<int> initScores(totalEdges);
    for (int e = 0; e < totalEdges; ++e) {
        initScores[e] = (int)B_edge[e];
    }
    K.initialize(initScores);

    vector<char> peeled(totalEdges, 0);
    vector<int>  peelSet;
    vector<pair<pair<int,int>,int>> contribBuf;
    int peeledCount = 0;

    while (peeledCount < totalEdges) {
        peelSet.clear();
        int e;
        while ((e = K.pop()) != -1) {
            if (!peeled[e]) {
                peeled[e] = 1;
                peelSet.push_back(e);
            }
        }
        peeledCount += (int)peelSet.size();

        updateE(neighbors, edgesList, edgeIdOf, incidentEdges,
                B_edge, peelSet, peelType, contribBuf);

        // #pragma omp parallel for schedule(dynamic)
        for (int idx = 0; idx < (int)peelSet.size(); ++idx) {
            int e1 = peelSet[idx];
            int u1 = edgesList[e1].first;
            int v1 = edgesList[e1].second;
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
    std::string line;

    std::vector<int> offsetsV, edgesV;
    std::vector<int> offsetsU, edgesU;

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
    for (int u = 0; u < v1; ++u) {
      int start = offsetsV[u], end = (u+1<v1 ? offsetsV[u+1] : e1);
      for (int j = start; j < end; ++j) {
        int w = v1 + edgesV[j];
        edgesList.emplace_back(u, w);
      }
    }

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
        scaleFactor = 1.0;
    }

    vector<vector<int>> neighbors(totalVertices);
    t1.start();
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
      edgeIdOf[{v,u}] = eid;
      incidentEdges[u].push_back(eid);
      incidentEdges[v].push_back(eid);
    }

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

    t1.start();
    vector<int> immediateNeighbors(totalVertices + 1);
    immediateNeighbors[0] = 0;
    for (int i = 0; i < totalVertices; ++i) {
        immediateNeighbors[i + 1] =
            immediateNeighbors[i] + (int)modifiedNeighbors[i].size();
    }
    int M = immediateNeighbors[totalVertices];

    vector<int> secondNeighbors(M);
    // #pragma omp parallel for schedule(dynamic)
    for (int u1 = 0; u1 < totalVertices; ++u1) {
        int modDeg = (int)modifiedNeighbors[u1].size();
        for (int i = 0; i < modDeg; ++i) {
            int v = modifiedNeighbors[u1][i];
            int cnt2 = 0;
            for (int w : neighbors[v]) {
                if (ranks[w] > ranks[u1]) {
                    ++cnt2;
                }
            }
            secondNeighbors[immediateNeighbors[u1] + i] = cnt2;
        }
    }

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
    vector<Wedge> W(totalWedges);
    t1.stop();
    t1.reportTotal("Initialising Wedge array");

    t1.start();
    // #pragma omp parallel for schedule(dynamic)
    for (int u1 = 0; u1 < totalVertices; ++u1) {
        int modDeg = modifiedNeighbors[u1].size();
        for (int i = 0; i < modDeg; ++i) {
            int v = modifiedNeighbors[u1][i];
            int base = (int)secondNeighborsPrefixSum[
                          immediateNeighbors[u1] + i];
            int cnt = 0;
            for (int w : neighbors[v]) {
                if (ranks[w] > ranks[u1])
                    W[base + cnt++] = {u1, w, v};
            }
        }
    }
    t1.stop();
    t1.reportTotal("Wedge retrieval");

    if (perType == "VERT") {
        vector<pair<pair<int,int>,int>> R;
        vector<int>                     F;
        if      (countType == "SORT")  aggregateWedgesSort(W, R, F);
        else if (countType == "HASH")  aggregateWedgesHash(W, R, F);
        else if (countType == "BATCH") aggregateWedgesBatch(W, R, F);
        else cmd.badArgument();

        vector<long long> B(totalVertices, 0LL);
        long long         totalButterflies = 0;

        // #pragma omp parallel for reduction(+ : totalButterflies)
        for (int i = 0; i < (int)R.size(); ++i) {
            int u1 = R[i].first.first;
            int u2 = R[i].first.second;
            int d  = R[i].second;
            long long c = (long long)d * (d - 1) / 2;
            // #pragma omp atomic
            B[u1] += c;
            // #pragma omp atomic
            B[u2] += c;
            totalButterflies += c;
        }

        // #pragma omp parallel for
        for (int i = 0; i < (int)R.size(); ++i) {
            int start = F[i], end = F[i+1];
            int d     = R[i].second;
            long long c = d - 1;
            for (int j = start; j < end; ++j) {
                int v = W[j].center;
                // #pragma omp atomic
                B[v] += c;
            }
        }
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
    else if (perType == "EDGE") {
        vector<pair<pair<int,int>,int>> R;
        vector<int>                     F;
        if      (countType == "SORT")  aggregateWedgesSort(W, R, F);
        else if (countType == "HASH")  aggregateWedgesHash(W, R, F);
        else if (countType == "BATCH") aggregateWedgesBatch(W, R, F);
        else cmd.badArgument();

        int T = 1; // omp_get_max_threads();
        vector<vector<pair<pair<int,int>,int>>> local_B_contribs(T);

        size_t totalContribs = 0;
        for (int i = 0; i < (int)R.size(); ++i) {
            totalContribs += 2 * (F[i+1] - F[i]);
        }

        for (int t = 0; t < T; ++t) {
            local_B_contribs[t].reserve(totalContribs / T + 1);
        }

        // #pragma omp parallel
        {
            int tid = 0; // omp_get_thread_num();
            auto& my_contribs = local_B_contribs[tid];

            // #pragma omp for schedule(dynamic)
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

        vector<pair<pair<int,int>,int>> B_contribs;
        B_contribs.reserve(totalContribs);
        for (int t = 0; t < T; ++t) {
            B_contribs.insert(B_contribs.end(),
                            local_B_contribs[t].begin(),
                            local_B_contribs[t].end());
        }

        vector<pair<pair<int,int>, long long>> edgeCounts;
        if      (countType == "SORT")  edgeCounts = aggregateKeyValueSort(B_contribs);
        else if (countType == "HASH")  edgeCounts = aggregateKeyValueHash(B_contribs);
        else if (countType == "BATCH") edgeCounts = aggregateKeyValueBatch(B_contribs);

        long long totalButterflies = 0;
        for (auto &ec : edgeCounts) {
            totalButterflies += ec.second;
        }
        cout << "Total number of butterflies = " << (long long)((totalButterflies / 4) * scaleFactor) << "\n";

        if (peelType != "NONE") {
            t1.start();

            vector<long long> B_edge(edgesList.size(), 0);
            for (const auto& ec : edgeCounts) {
                int eid;
                auto it = edgeIdOf.find(ec.first);
                if (it != edgeIdOf.end()) {
                    eid = it->second;
                }
                else {
                    auto it2 = edgeIdOf.find({ec.first.second, ec.first.first});
                    if (it2 != edgeIdOf.end())
                        eid = it2->second;
                    else continue;
                }
                B_edge[eid] = ec.second;
            }

            peelE(
                (int)edgesList.size(),
                neighbors,
                edgesList,
                edgeIdOf,
                incidentEdges,
                B_edge,
                peelType
            );

            t1.stop();
            t1.reportTotal("Time To Peel Edges");
        }

        startTime.stop();
        startTime.reportTotal("Total time to run");
    }

    return 0;
}