# PDC_Project

Repository for Parallel and Distributed Computing Project where the research paper https://epubs.siam.org/doi/pdf/10.1137/1.9781611976021.2 will be analyzed and implemented for better understanding how parallel and distributed computers function and how they can be used to achieve significant speedup when processing large datasets.

## Introduction

In bipartite graphs, a Butterfly is a four-vertex subgraph with two vertices on each side. It represents the following edges:

U1 ---- V1
U1 ---- V2
U2 ---- V1
U2 ---- V2

When this is plotted properly this will produce the shape of a Butterfly, however some people may also refer to it as an Hourglass. Butterfly computations are fundamental for detecting patterns and dense local structures, showing overlap, similarity, and community behaviour in a network. This is further used for ranking algorithms, visualizations, and compressions of graphs. They provide a rich structural insight into how entities in a bipartite system relate and interact. A list of real world examples is as follows:

- Movie platforms: Users who co-watch the same films
- E-commerce: Products bought together by groups of users
- Academic networks: Co-authorship patterns in publications
- Social media: Common interactions (hashtags, follows, posts)

This specific paper deals how Bipartite Graphs and Butterflies present within them can be properly analyized and processed through the use of a Parallel Framework known as ParButterfly, which is meant to make Butterfly computations scalable and efficient across multi-core systems. It addresses several key tasks:

- Global Butterfly Counting
- Per-Vertex/Per-Edge Counting
- Peeling Algorithms

Instead of using a one-size-fits-all approach, several parallel strategies are tailored to each task, such as sorting-based wedge grouping, hashing and batching, degree-aware vertex ordering, and graph sparsification for approximatinos.

The results of this were stated in the paper as follows:
- 6.3x to 13.6x speedups over best sequential implementations
- 349.6x to 5169x speedup over best parallel implementations
- 7.1x to 38.5x self-relative speedups
- Up to 1.7x additional speedup using cache-optimization

This repository will be used to implement a custom version of this solution using OpenMP, OpenCL, MPI, and METIS.
