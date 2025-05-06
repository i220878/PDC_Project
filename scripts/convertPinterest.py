#!/usr/bin/env python3
import sys

def main(field_path, edges_path):
    # 1) Read field.dat → id2type ('u' or 'i')
    id2type = {}
    max_id = -1
    with open(field_path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            tok = line.split()
            vid = int(tok[0])
            typ = tok[1]
            id2type[vid] = typ
            if vid > max_id: max_id = vid

    # 2) Assign new indices for users (V‐side) and items (U‐side)
    id2v = {}  # old userID -> [0..nv-1]
    id2u = {}  # old itemID -> [0..nu-1]
    nv = nu = 0
    for vid in range(max_id+1):
        t = id2type.get(vid)
        if t == 'u':
            id2v[vid] = nv
            nv += 1
        elif t == 'i':
            id2u[vid] = nu
            nu += 1

    # 3) Read edges.dat → build bipartite adjacency lists
    adjV = [[] for _ in range(nv)]  # V‐side → list of U‐indices
    adjU = [[] for _ in range(nu)]  # U‐side → list of V‐indices
    with open(edges_path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split()
            if len(parts) < 3:
                continue
            a, b, _ = parts
            pa, xa = a[0], int(a[1:])
            pb, xb = b[0], int(b[1:])
            if pa == 'u' and pb == 'i':
                u0, v0 = xa, xb
            elif pa == 'i' and pb == 'u':
                u0, v0 = xb, xa
            else:
                continue
            if u0 not in id2v or v0 not in id2u:
                continue
            u = id2v[u0]
            v = id2u[v0]
            adjV[u].append(v)
            adjU[v].append(u)


    # 4) Deduplicate & sort each list
    for L in adjV:
        L[:] = sorted(set(L))
    for L in adjU:
        L[:] = sorted(set(L))

    # 5) Build CSR‐style offsets & edge arrays
    offV = [0]
    for L in adjV:
        offV.append(offV[-1] + len(L))
    edgeV = [v for L in adjV for v in L]

    offU = [0]
    for L in adjU:
        offU.append(offU[-1] + len(L))
    edgeU = [u for L in adjU for u in L]

    mv = offV[-1]
    mu = offU[-1]

    # 6) Emit AdjacencyHypergraph
    out = sys.stdout
    out.write("AdjacencyHypergraph\n")
    out.write(f"{nv}\n{mv}\n{nu}\n{mu}\n")
    for x in offV[:-1]: out.write(f"{x}\n")
    for v in edgeV:    out.write(f"{v}\n")
    for x in offU[:-1]: out.write(f"{x}\n")
    for u in edgeU:    out.write(f"{u}\n")

if __name__ == "__main__":
    if len(sys.argv)!=3:
        print("Usage: to_hypergraph.py field.dat edges.dat", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
