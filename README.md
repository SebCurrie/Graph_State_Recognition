# Graph_Code
Code written to analyse whether a given input state is a qubit/qudit graph/hypergraph.


The code is split into methods for qubits and methods for qudits. Within each of these, there are four files.
The graph_finder and graph_methods files perform a brute force stabiliser check on the state vector to work out if a given input state is a graph state. If it is graph_methods will return the graph. These methods are inefficient and so we implement a nicer algorithm.

The hypergraph_finder files will use this more efficient algorithm and hypergraph_methods will return the graph. The algorithm for qubits is given in https://arxiv.org/abs/1211.5554. We generalised the algorithm to extend to hypergraphs of qudits. The theory of this work will be posted here soon.

Full README to follow. Several minor issues still to be fixed and some functionality to be added.


