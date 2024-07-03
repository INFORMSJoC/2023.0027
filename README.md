[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# CacheTest

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The data in this repository are a snapshot of the data
that were used in the research reported on in the paper 
[This is a Template](https://doi.org/10.1287/ijoc.2019.0000) by Zhong et. al.. 
The snapshot is based on 
[this SHA](https://github.com/tkralphs/JoCTemplate/commit/f7f30c63adbcb0811e5a133e1def696b74f3ba15) 
in the development repository. 

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2019.0000

https://doi.org/10.1287/ijoc.2019.0000.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{CacheTest,
  author =        {H. Zhong, F. Mahdavi Pajouh, S. Butenko, O.A. Prokopyev},
  publisher =     {INFORMS Journal on Computing},
  title =         {{On interdicting dense clusters in a network}},
  year =          {2024},
  doi =           {10.1287/ijoc.2019.0000.cd},
  url =           {https://github.com/INFORMSJoC/2019.0000},
  note =          {Available for download at https://github.com/INFORMSJoC/2019.0000},
}  
```

## Description

Data set used in “On interdicting dense clusters in a network” paper. It contains 20 randomly-generated networks and 38 real-life networks. Each file contains 1+X+Y lines, where X and Y are the number of vertices and the number of edges in the network, respectively. It starts with two numbers indicating the number of vertices X and the number of edges Y. Then, it follows with X lines of two numbers. Note that we do not assign indices to the vertices. The index of a vertex is the same as the line number it appears on. The first number in each line is the weight of a vertex and the second number is its blocking cost. The last Y lines contain three numbers. The first and second numbers in each line are the head and the tail of an edge. Note that the edges are undirected. The third number is the blocking cost of an edge.
