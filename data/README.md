# Data for replication

Folder "random" contains 20 randomly-generated networks and folder "real life" contains 38 real-life networks. Each file contains 1+X+Y lines, where X and Y are the number of vertices and the number of edges in the network, respectively. It starts with two numbers indicating the number of vertices X and the number of edges Y. Then, it follows with X lines of two numbers. Note that we do not assign indices to the vertices. The index of a vertex is the same as the line number it appears on. The first number in each line is the weight of a vertex and the second number is its blocking cost. The last Y lines contain three numbers. The first and second numbers in each line are the head and the tail of an edge. Note that the edges are undirected. The third number is the blocking cost of an edge.

For real-life instances, we chose 38 graphs from different categories of Network Repository cite below.

R. A. Rossi and N. K. Ahmed. The network data repository with interactive graph analytics and visualization. In AAAI,
2015. URL https://networkrepository.com.
