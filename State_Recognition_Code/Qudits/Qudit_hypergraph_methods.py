import networkx as nx
import matplotlib.pyplot as plt
import hypernetx as hnx
import copy
import math



#This file contains functions for going between qudit hypergraphs and stabilisers
#The input to the get_hypergraph_from_edges() function takes a list of edges as input




#Function which prints the hypernetx graph
def hypernetx_graph_printer(outputG,stabilisers):
    G = outputG
    
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(111)
    ax.set_title('Graph - Shapes', fontsize=10)
    hnx.drawing.rubber_band.draw(G, pos=None, with_color=True, with_node_counts=False, with_edge_counts=False, layout_kwargs={}, ax=None, edges_kwargs={}, nodes_kwargs={}, edge_labels_kwargs={}, node_labels_kwargs={}, with_edge_labels=True, with_node_labels=True, label_alpha=0.35)
    plt.tight_layout()
    plt.xlabel('HyperNetX output graph')
    #plt.savefig("Graph.pdf", format="PDF")
    
    return



#### USING HYPERNETWORKX NOT NETWORKX
#This method will not properly print any nodes which are not entangled!!
#This function takes some edges, turns them into the correct format for hypernetx to print
# and then passes it to the printer function to display
def get_hypergraph_from_edges(edges,n):
    entities={}
    
    #This is all just to get it to display nodes with no links
    #Doesnt even work properly, displays them as self loops rather then nodes. need to tidy up
    for i in range(n):
        linked = False

        if len(edges)==0:
            entities[str(i)]=[]
        else:    
            for j in range(len(edges)):
                if i in edges[j]:
                    linked=True                 
            if linked==False:    
                entities[str(i)]=[]
                
    
    #Convert edges into appropriate dictionary format for hypernetx
    for i in range(len(edges)):
         entities[str(i)]=edges[i]
         

    output_graph = hnx.Hypergraph(entities)
    get_hypergraph_stabilisers(edges,n)
    hypernetx_graph_printer(output_graph,edges)

    return output_graph



#Returns a list of hypergraph stabilisers given the hypergraph edges
def get_hypergraph_stabilisers(edges,n):
    stab_list=[['X'+str([i])] for i in range(n)]
    
    for i in range(n):
       
        for j in range(len(edges)):
            
            if i in edges[j]:
                temp=[]
                temp.extend(copy.deepcopy(edges[j]))
                
                for x in range(len(temp)):
                    if temp[x]==i:
                        del(temp[x])
                        break
                
                stab_list[i].append('CZ'+str(temp))
            else:
                continue
                
                

    print('\n The stabilisers of this hypergraph are: \n',stab_list,'\n')

    return stab_list


