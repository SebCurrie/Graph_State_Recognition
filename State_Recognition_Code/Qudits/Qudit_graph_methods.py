import networkx as nx
import matplotlib.pyplot as plt
import math

#Methods for plotting graph states of qudits



#NETWORKX METHODS ARE BELOW THIS ====================================================================
#Unfortunately Graph State compass uses networkx which has no hypergraph capability.
#The following methods will work for graphs of qudits, no hypergraphs of qudits.


#Function taken from Graph-State-Compass
def is_prime(a):
    if a < 2:
        return False
    for x in range(2, int(math.sqrt(a)) + 1):
        if a % x == 0:
            return False
    return True


#Function taken from Graph-State-Compass
def nx_create_prime_graph(w_edges, prime):
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

#Function taken from Graph-State-Compass
def nx_create_prime_power_graph(w_edges, prime, power):
    """ Creates a weighted graph representing a prime qudit graph state """
    nx_wg = nx_create_prime_graph(w_edges, prime)
    nx_wg.power, nx_wg.dimension = power, int(prime ** power)
    fam_labels = list(set([n for n, i in nx_wg.nodes()]))
    print(fam_labels)
    nx_wg.families = len(fam_labels)
    fam_nodes = [(n, i) for n in fam_labels for i in range(power)]
    # Adds any nodes that weren't in the edge list
    nx_wg.add_nodes_from(fam_nodes)
    
    return nx_wg
    
    
    

def nx_get_qudit_graph_from_stabilisers(stab_list,d):
    edge=[]
    edges=[]
    nodes=[]
    
    for i in range(len(stab_list)):
        x_node=0
        z_node=[]
        weight=[]
        
        nodes.append(i)
        
        for j in range(len(stab_list)):
            temp=stab_list[i][j]
            
            if temp=='I':
                continue
            elif temp[0]=='X':
                x_node=j
                
                
            elif temp[0]=='Z':
                z_node.append(j)
                weight.append(int(temp[2]))
        
        for a in range(len(z_node)):       
             edge=(x_node,z_node[a],weight[a])
             edges.append(edge)
    
    graph=nx_create_prime_graph(edges,d)
    nx_qudit_graph_printer(graph,d)
    
    return graph
    
    
    
def nx_qudit_graph_printer(outputG,dimension):
    G = outputG
    
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(111)
    ax.set_title('Qudit Graph states. Dimension='+str(dimension), fontsize=20)
    plt
    labels = nx.get_edge_attributes(G,'weight')
    
    pos = nx.spring_layout(G)
    
    nx.draw_networkx_nodes(G, pos, node_size=700)
    nx.draw_networkx_edges(G, pos, edgelist=G.edges, width=6)
    nx.draw_networkx_edge_labels(G,pos,edge_labels=labels, font_size=20)
    nx.draw_networkx_labels(G, pos, font_size=20, font_family='sans-serif')
 
    plt.tight_layout()
    plt.xlabel('Qudit Graph State')
    #plt.savefig("Graph.pdf", format="PDF")
    
    
    return
