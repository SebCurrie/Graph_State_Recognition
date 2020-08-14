import networkx as nx
import matplotlib.pyplot as plt
from IPython.display import display
import numpy as np
from halp import undirected_hypergraph as hp
import hypernetx as hnx


def hypernetx_graph_printer(outputG,stabilisers):
    G = outputG
    
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(111)
    ax.set_title('Graph - Shapes', fontsize=10)

    

    hnx.drawing.rubber_band.draw(G, pos=None, with_color=True, with_node_counts=False, with_edge_counts=False, layout_kwargs={}, ax=None, edges_kwargs={}, nodes_kwargs={}, edge_labels_kwargs={}, node_labels_kwargs={}, with_edge_labels=True, with_node_labels=True, label_alpha=0.35)
    plt.tight_layout()
    plt.xlabel('sd')
    plt.savefig("Graph.pdf", format="PDF")


    return



#### USING HALP NOT NETWORKX
def get_hypergraph_from_edges(edges,n):
    nodes=[]
    entities={}
    for i in range(len(edges)):
         entities[str(i)]=edges[i]
         

    output_graph = hnx.Hypergraph(entities)
    hypernetx_graph_printer(output_graph,edges)
get_hypergraph_from_edges([[1,2],[2,3],[1,3],[1,2,3],[1]],3)






#USING NETWORKX------No hypergraph functionality

def get_graph_from_stabilisers(input_stabs):
    edge=[]
    edges=[]
    nodes=[]
    
    for i in range(len(input_stabs)):
        x_node=0
        z_node=[]
        
        
        nodes.append(i)
        
        for j in range(len(input_stabs)):
            
            if input_stabs[i][j]=='X':
                x_node=j
                
            if input_stabs[i][j]=='Z':
                z_node.append(j)
                
        for a in z_node:       
             edge=(x_node,a)
             edges.append(edge)
        
   
    
    edges_checked=edge_redundancy(edges)
    
    output_graph = nx.Graph()
    output_graph.add_nodes_from(nodes)
    output_graph.add_edges_from(edges_checked)
    
    graph_printer(output_graph,input_stabs)




def get_stabilisers_from_graph(graph_obj):
    
    no_nodes=len(graph_obj.nodes)
    stabilisers = [[] for x in range(no_nodes)]
    
  
    
    for i in graph_obj.nodes:
        stab_temp=[[] for x in range(no_nodes)]
        stab_temp[i]='X'
        
        for j in range(no_nodes):
           if j==i:
               continue
           elif graph_obj.has_edge(i, j):
               stab_temp[j]='Z'
           else:
               stab_temp[j]='I'
        
        stabilisers[i]=stab_temp
        graph_printer(graph_obj,stabilisers)
    return stabilisers
 
           
    


def edge_redundancy(edge_list):

    for i,j in edge_list:
            if edge_list[i][0]==edge_list[j][1]:
                del edge_list[j]
    
    return edge_list




    

def graph_printer(outputG,stabilisers):
    G = outputG
    
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(111)
    ax.set_title('Graph - Shapes', fontsize=10)
    
    # some math labels
    labels={}
    for i in G.nodes:
        labels[i]=str(i)
        fig.text((0.1+i*0.2), 0.08, stabilisers[i], ha='center', fontsize=12)

    pos = nx.spring_layout(G)
    nx.draw_networkx_labels(G,pos,labels,font_size=16)
    nx.draw(G, pos, node_size=1500, node_color='red', font_size=8, font_weight='bold')
    
    
    plt.tight_layout()
    plt.xlabel('sd')
    plt.savefig("Graph.pdf", format="PDF")

  
    
# Create the 4 qubit line graph from Jeremys Paper
#Note that the stabilisers are different under relableling of the nodes
# =============================================================================
# input_edges = [(2, 0), (0, 1), (1, 3)]
# inputG = nx.Graph()
# inputG.add_edges_from(input_edges)
# 
# =============================================================================

#4 qubit star from Jeremeys Paper
#input_edges = [(2, 3), (0, 3), (1, 3)]
#inputG = nx.Graph()
#inputG.add_edges_from(input_edges)

#stab_list=get_stabilisers(inputG)
#stab_list=[['X', 'I', 'I'], ['I', 'X', 'I'], ['I', 'I', 'X']]
    
    
#stab_list=[['X', 'Z', 'Z', 'I'], ['Z', 'X', 'I', 'Z'], ['Z', 'I', 'X', 'I'], ['I', 'Z', 'I', 'X']]
#get_graph(stab_list)
