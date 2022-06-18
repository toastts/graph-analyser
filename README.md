# Graphs
implementations of graphs and functions using the properties of graphs

## Prerequisites

In order to run the Markov Clustering algorithm the Eigen 3 Library, in which version 3.3.7 is the latest stable version. This can be downloaded from [Eigen Main Page](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is the only necessary library that must be downloaded in order to run code.

## Data

All the data used were various types of graphs such as directed and weighted, undirected and unweighted, etc. Imported all data from [SNAP Database](https://snap.stanford.edu/data/). The methods interepret interesting properties of graphs such as being cyclic or not, clustering using multiple methods, indegree and outdegree, and max flow. 

## Code

code is made up of various hpp files:

#### [edge.hpp](https://github.com/BASIS-Bradley/TeamJunior/blob/master/edge.hpp)
  * Makes a single edge object that contains source, destination, and weight info with a comparator that allows for minheaps based on edge weight to be created in analyser methods

#### [node.hpp](https://github.com/BASIS-Bradley/TeamJunior/blob/master/node.hpp)
  * Most basic object for each point in the graph containing the number of the node
  
#### [point.hpp](https://github.com/BASIS-Bradley/TeamJunior/blob/master/point.hpp)
  * Object representing each point of the graph that creates a node as well as retrieves it and its weight
  
#### [graph.hpp](https://github.com/BASIS-Bradley/TeamJunior/blob/master/graph.hpp)
  * Object that creates and prints a collection of the points in the graph interrelated by edges as well as retrieves its properties
  
#### [analyser.hpp](https://github.com/BASIS-Bradley/TeamJunior/blob/master/analyser.hpp)
  * Contains all but one of the analysis funtions that can be run on the graphs, which writes the results to text files 
  
#### [mcl.hpp](https://github.com/BASIS-Bradley/TeamJunior/blob/master/mcl.hpp)
  * Contains Markov Clustering Algorithm
  
#### [main.hpp](https://github.com/BASIS-Bradley/TeamJunior/blob/master/main.cpp)
  * Contains the method to create and import graphs as well as run funtions on those graphs

## Running the Code

Create an executable of [main.hpp file](https://github.com/BASIS-Bradley/TeamJunior/blob/master/main.cpp) by typing and then run that executable. This will write the results of the analysis functions on the graphs into multiple text files.
