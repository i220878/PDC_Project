#!/usr/bin/env python3
"""
Convert an undirected edge‐list with comments (DBLP coauthor graph)
into the “AdjacencyHypergraph” bipartite‐CSR format:

  Part A: treat each author as a V‐node
  Part B: treat each original edge as a U‐node (incidence)
"""
import sys

def usage():
    print("Usage: dblp_to_csr.py INFILE OUTFILE")
    sys.exit(1)

if len(sys.argv) != 3:
    usage()
INFILE  = sys.argv[1]
OUTFILE = sys.argv[2]

# 1) Read edge‐list, skip comments, track max node
edges = []
max_id = -1
with open(INFILE) as f:
    for line in f:
        line = line.strip()
        if not line or line[0] == '#':
            continue
        u, v = map(int, line.split())
        edges.append((u, v))
        max_id = max(max_id, u, v)

V = max_id + 1          # number of authors = |V|
E = len(edges)          # number of original edges = |U|

print(f"Read {E} edges over {V} authors")

# 2) Build degree of each V‐node (how many incident edges)
degV = [0] * V
for (u, v) in edges:
    degV[u] += 1
    degV[v] += 1

# 3) Build V→U CSR (offsetsV, edgesV)
offsetsV = [0]*V
for i in range(1, V):
    offsetsV[i] = offsetsV[i-1] + degV[i-1]
# total entries = sum(degV) = 2*E
e1 = 2*E
edgesV = [0]*e1
# fill in neighbors
cur = offsetsV.copy()
for eid, (u, v) in enumerate(edges):
    edgesV[cur[u]] = eid
    cur[u] += 1
    edgesV[cur[v]] = eid
    cur[v] += 1

# 4) Build U→V CSR (offsetsU, edgesU)
# each U‐node (original edge) has exactly two endpoints
offsetsU = [2*i for i in range(E)]
e2 = 2*E
edgesU   = [0]*(e2)
for eid, (u, v) in enumerate(edges):
    edgesU[2*eid]   = u
    edgesU[2*eid+1] = v

# 5) Write out in AdjacencyHypergraph format
with open(OUTFILE, 'w') as o:
    o.write("AdjacencyHypergraph\n")
    o.write(f"{V}\n")    # nv
    o.write(f"{e1}\n")   # eV = |edgesV|
    o.write(f"{E}\n")    # nu
    o.write(f"{e2}\n")   # eU = |edgesU|
    # offsetsV
    for x in offsetsV:
        o.write(f"{x}\n")
    # edgesV
    for x in edgesV:
        o.write(f"{x}\n")
    # offsetsU
    for x in offsetsU:
        o.write(f"{x}\n")
    # edgesU
    for x in edgesU:
        o.write(f"{x}\n")

print(f"Wrote CSR bipartite to {OUTFILE}:")
print(f"  nv={V}, eV={e1}, nu={E}, eU={e2}")
