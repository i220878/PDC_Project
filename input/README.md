Input Format for ParButterfly
-----------
The input format of bipartite graphs is based on the adjacenecy graph format
from the Problem Based Benchmark Suite and Ligra. The adjacency graph format
starts with a sequence of offsets one for each vertex in one bipartition V,
 followed by a sequence of directed edges ordered by their source
 vertex. The offset for a vertex i refers to the location of the start
 of a contiguous block of out edges for vertex i in the sequence of
 edges. The block continues until the offset of the next vertex, or
 the end if i is the last vertex. All vertices and offsets are 0 based
 and represented in decimal. This then repeats for the second bipartition U.

 Let nv and nu denote the number of vertices in bipartitions V and U respectively. 
 Let mv and mu denote the number of edges (mv = mu). 
 The specific format is as follows:

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