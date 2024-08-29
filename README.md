[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# On interdicting dense clusters in a network

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The data in this repository were used in the research reported on in the paper 
[On interdicting dense clusters in a network](https://doi.org/10.1287/ijoc.2023.0027) by Zhong et al. 

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0027

https://doi.org/10.1287/ijoc.2023.0027.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{CacheTest,
  author =        {Zhong, Haonan and Mahdavi Pajouh, Foad and Butenko, Sergiy and Prokopyev, Oleg A.},
  publisher =     {INFORMS Journal on Computing},
  title =         {{On interdicting dense clusters in a network}},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0027.cd},
  url =           {https://github.com/INFORMSJoC/2023.0027},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0027},
}  
```

## Description

Source codes and data set used in “On interdicting dense clusters in a network” paper.  Data set contains 20 randomly-generated networks and 38 real-life networks. Each file contains 1+X+Y lines, where X and Y are the number of vertices and the number of edges in the network, respectively. It starts with two numbers indicating the number of vertices X and the number of edges Y. Then, it follows with X lines of two numbers. Note that we do not assign indices to the vertices. The index of a vertex is the same as the line number it appears on. The first number in each line is the weight of a vertex and the second number is its blocking cost. The last Y lines contain three numbers. The first and second numbers in each line are the head and the tail of an edge. Note that the edges are undirected. The third number is the blocking cost of an edge.
