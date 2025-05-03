# read your V-side offsets and edges into Python lists
with open("/home/harth/Desktop/Work/School/PDC/Project/parbutterfly/inputs/rMatGraph_J_5_100.txt") as f:
    lines = f.read().split()
# skip header magic and the four ints
it = iter(lines)
assert next(it) == "AdjacencyHypergraph"
nv, mv, nu, mu = map(int, (next(it), next(it), next(it), next(it)))
# read V-offsets
offsetv = [int(next(it)) for _ in range(nv)]
# read V→U edges
edgev   = [int(next(it)) for _ in range(mv)]

# build U adjacency lists
Uadj = [[] for _ in range(nu)]
for v in range(nv):
    start = offsetv[v]
    end   = (offsetv[v+1] if v+1<nv else mv)
    for ei in range(start, end):
        u = edgev[ei]
        Uadj[u].append(v)

# now build CSR U-side
offsetu = [0]*nu
for u in range(1,nu):
    offsetu[u] = offsetu[u-1] + len(Uadj[u-1])

edgeu = []
for u in range(nu):
    edgeu.extend(Uadj[u])

assert len(edgeu) == mu == mv

# print out the complete file
print("AdjacencyHypergraph")
print(nv);  print(mv)
print(nu);  print(mu)
for x in offsetv: print(x)
for x in edgev:   print(x)
for x in offsetu: print(x)
for x in edgeu:   print(x)
