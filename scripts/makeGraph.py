#!/usr/bin/env python3
import sys, argparse, random

def main():
    p = argparse.ArgumentParser(
        description="Generate a random bipartite graph in AdjacencyHypergraph format"
    )
    p.add_argument('nv', type=int,
                   help='number of vertices in the V partition')
    p.add_argument('nu', type=int,
                   help='number of vertices in the U partition')
    p.add_argument('m', type=int,
                   help='number of bipartite edges to pick at random')
    p.add_argument('--seed', type=int, default=None,
                   help='random seed (default: None)')
    args = p.parse_args()

    nv, nu, m = args.nv, args.nu, args.m
    if args.seed is not None:
        random.seed(args.seed)

    max_edges = nv * nu
    if m < 0 or m > max_edges:
        sys.exit(f"error: m must be between 0 and {max_edges}")

    # 1) sample m distinct edges (v,u) with v in [0..nv), u in [0..nu)
    edges = set()
    while len(edges) < m:
        v = random.randrange(nv)
        u = random.randrange(nu)
        edges.add((v,u))

    # 2) build adjacency lists V->U and U->V
    V_adj = [[] for _ in range(nv)]
    U_adj = [[] for _ in range(nu)]
    for (v,u) in edges:
        V_adj[v].append(u)
        U_adj[u].append(v)

    # 3) CSR for V side
    offsetsV = [0]*nv
    edgesV   = []
    for i in range(nv):
        offsetsV[i] = len(edgesV)
        # append all incident U-indices
        for u in V_adj[i]:
            edgesV.append(u)
    mv = len(edgesV)  # should equal m

    # 4) CSR for U side
    offsetsU = [0]*nu
    edgesU   = []
    for j in range(nu):
        offsetsU[j] = len(edgesU)
        for v in U_adj[j]:
            edgesU.append(v)
    mu = len(edgesU)  # should also equal m

    # 5) print in AdjacencyHypergraph format
    out = sys.stdout
    out.write("AdjacencyHypergraph\n")
    out.write(f"{nv}\n")
    out.write(f"{mv}\n")
    out.write(f"{nu}\n")
    out.write(f"{mu}\n")

    # V-offsets
    for off in offsetsV:
        out.write(f"{off}\n")
    # V-edges
    for e in edgesV:
        out.write(f"{e}\n")

    # U-offsets
    for off in offsetsU:
        out.write(f"{off}\n")
    # U-edges
    for e in edgesU:
        out.write(f"{e}\n")

if __name__ == "__main__":
    main()
