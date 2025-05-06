# Implementation

Currently, sequential and parallel implementations have been done:
1. Degree-based Ranking
2. Sorting for Wedge aggregation
3. SEQUENTIAL butterfly counting, NO atomic adds or batching yet.
4. NO Peeling / decomposition implemented.
5. NO Approximation implemented.

Simply compile the file with `make` and run `./SequentialBipartiteButterflies` for the sequential version, or `./ParallelBipartiteButterflies` for the Bipartite Implementation.
Further parallel optimisations should be made

The input format is as follows:

AdjacencyHypergraph  
&lt;nv>  
&lt;mv>  
&lt;nu>  
&lt;mu>  
&lt;offsetv(0)>   
&lt;offsetv(1)>  
...  
&lt;offsetv(nv-1)>  
&lt;edgev(0)>  
&lt;edgev(1)>  
...  
&lt;edgev(mv-1)>   
&lt;offsetu(0)>  
&lt;offsetu(1)>  
...  
&lt;offsetu(nv-1)>  
&lt;edgeu(0)>  
&lt;edgeu(1)>  
...  
&lt;edgeu(mv-1)>

For more details, check the README.md in parbutterfly/

# Running the code

The following input format is expected to run the code (and is inspired by the paper):

The application ParallelBipartiteButterflies takes a graph as input as shown in the previous section.
It also takes parameters for counting and peeling, as follows:

| Parameter    | Default     | Description                                      |
| ---------    | -------     | ------------------------------------------------ |
| `-countType` | BATCH       | (SORT, HASH, BATCH) Wedge/butterfly aggregation type for counting |
| `-rankType`  | ADEG        | (SIDE, DEG, or ADEG) Vertex ordering . The prefix A means approximate.                |
| `-peelType`  | NONE        | (SORT, HASH, BATCH, NONE) Butterfly aggregation type for peeling. NONE means don't execute peeling. |
| `-per`       | VERT        | (VERT, EDGE) Denotes whether to count per vertex (VERT), or per edge (EDGE) |
| `-sparseType`| NONE        | Sparsification options (EDGE, COLOR, NONE) |
| `-d`         | 5           | If -sparseType is EDGE, then 1/d is the probability of keeping an edge. If -sparseType is COLOR, then d is the number of colors. |

# NOTE
The first 5 lines of the input file MUST be as shown above, with the first being "AdjacencyHypergraph",
and the following four being
nv (number of vertices for first bipartition)
mv (number of edges for first bipartition)
nu (number of vertices for second bipartition)
mu (number of edges for second bipartition)

Note that the provided file in the default repository of parbutterfly/ was not correct and had to be altered. Afterwards, this code and both the parbutterfly code gives the same output for `..input/full_graph.txt` at 712 butterflies, though the parbutterfly implementation is orders of magnitude faster.

To run the parbutterfly implementation, navigat to parbutterfly/apps/hyper and write the following command: `make OPENMP=1 -j`

followed by `./BipartiteButterflies -countType SORT -rankType DEG -peelType NONE ../../inputs/full_graph.txt` with the arguments altered if required to the options given in `parbutterfly/README.md`.