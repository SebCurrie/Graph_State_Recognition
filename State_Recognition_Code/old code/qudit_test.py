import itertools
import Graph_Printers as gp
import math
import networkx as nx
import matplotlib.pyplot as plt
from IPython.display import display
import numpy as np
import hypernetx as hnx
import copy


def is_prime(a):
    if a < 2:
        return False
    for x in range(2, int(math.sqrt(a)) + 1):
        if a % x == 0:
            return False
    return True


def create_prime_graph(w_edges, prime):
    """ Creates a weighted graph representing a prime qudit graph state """
    if not is_prime(prime):
        raise Exception("Graph state must be prime-dimensional")
    us, vs, ws = zip(*w_edges)
    if max(ws) >= prime or max(ws) < 0:
        raise Exception("Weights must be 0 <= w < p ")
    nx_wg = nx.Graph()
    nx_wg.add_weighted_edges_from(w_edges)
    nx_wg.prime = prime
    nx_wg.power = 1
    nx_wg.dimension = prime
    
    return nx_wg


def create_prime_power_graph(w_edges, prime, power):
    """ Creates a weighted graph representing a prime qudit graph state """
    nx_wg = create_prime_graph(w_edges, prime)
    nx_wg.power, nx_wg.dimension = power, int(prime ** power)
    fam_labels = list(set([n for n, i in nx_wg.nodes()]))
    print(fam_labels)
    nx_wg.families = len(fam_labels)
    fam_nodes = [(n, i) for n in fam_labels for i in range(power)]
    # Adds any nodes that weren't in the edge list
    nx_wg.add_nodes_from(fam_nodes)
    
    return nx_wg

#def get_qudit_state_vec_from_graph()

w_edges = [(0, 1, 1), (1, 2, 2)]                #(node1,node2,weight)
qudit_graph=create_prime_graph(w_edges,3)

w_edges = [(0, 1,1)]                #(node1,node2,weight)
qudit_graph=create_prime_power_graph(w_edges,2,2)
print()
gp.qudit_graph_printer(qudit_graph,4)





