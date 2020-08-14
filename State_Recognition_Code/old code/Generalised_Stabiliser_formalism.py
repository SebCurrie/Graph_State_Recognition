import itertools
import Graph_Printers as gp
import math
import networkx as nx
import matplotlib.pyplot as plt
from IPython.display import display
import numpy as np
import hypernetx as hnx
import copy
import cmath


#THIS CODE SHOULD NOT BE USED --- LOOK FOR SPARSE MATRIX IMPLEMENTATION IN QUDIT STABILSERS SPARSE

def generate_generalised_pauli_matrices(d):
    w=cmath.exp(2*math.pi*1j/d)
    
    #Generate the "sub" matrices before combining to create elements of the pauli group
    #Will return d matrices for b=0,1,2...d-1
    zs=[]
    for l in range(d):
        z=np.array(np.zeros((d,d)),dtype=complex)
        
        for t in range(d):
            z[t][t]=w**(l*t)
        zs.append(z)
        
    
    #genenrate the "sub" x matrices
    #These x's are always just identiy matrices with the rows swapped in a procedural way
    #Write one out to check this algorithm works
    #Will generate d matrices for a=0,1,2...d-1
    #Easier way to do this would be to just take powers the initial x matrix
    xs=[]
    I=np.identity(d)
    for h in range(d):
        
        x=np.array(np.zeros((d,d)))
        for f in range(d):
            x[f]=I[(f+(d-h))%(d)]
        xs.append(x)
        
    #generate ana array of the appropraite roots of unity w**(c)
    ws=[]
    for m in range(d):
        ws.append(w**m)
   
    return xs,zs,ws


#Function that checks all if each combination of generalised pauli is a stabiliser
#To get the minimal generating set of the stabiliser group it is sufficient to only test the
# elements which have X(1). The rest of the elements X(2)...X(d-1) are just powers of the X(1) stabiliser
#The X(1) elements also contain the correct edge weight information.
    
def qudit_brute_force_stabiliser_check(input_state,d,n):
    xs,zs,ws=generate_generalised_pauli_matrices(d)
    stabs=[]
    
     #This generates a cartesian product of all possible combs of Z and I
    #doing like this allows us to just input n. We then copy the list
    #multiple times and insert X into a different position on each list 
    elements=['I']
    for i in range(d):
        elements.append('Z('+str(i)+')')
        
    
    combo=itertools.product(elements, repeat=n-1)
   
    combos=[list(p) for p in combo]
   
    for i in range(n):
        stabs.extend(copy.deepcopy(combos))
    
   
    for i in range(n):
        for j in range(len(combos)):
           stabs[i*(len(combos))+j].insert(i,'X(a)')
           
    print('The set of all possible stabilsers to be tested for n=',n,'qudit dimension d=',d ,'is: \n', stabs, '\n')
    
    
    #This generates a symbolic list of all the stabilisers that need to be checked
    #     #Now we are in the qudit case there are multiple X matrices, and multiple Z matrices to check
#     #and we must check all possible combs
#     # i.e X(0)Z(0), X(0)Z(1), X(1),Z(0), X(1)Z(1)
#     #     Z(0)X(0), Z(1)X(0),  Z(0)X(1)  , Z(1)X(1)
    qudit_stabs=[]
    for i in range(len(stabs)):
        temp=stabs[i]
        for k in range(n):
            if temp[k]=='X(a)':
                qudit_stab='X('+str(1)+')'
                temp[k]=qudit_stab
                qudit_stabs.append(copy.deepcopy(temp))
    

   
    print(' As a,b go from 0 to',d-1, 'this is the set: \n' ,qudit_stabs)
    
    
    
#     # Now given the list of all possible combinations we need to test, we generate a matrix for each
#     #combination and multiply it by the input state vector. If its an eigenvector we append it to the list of 
#     # checked stabilisers
    stab_mats=[]
    stabiliser_list=[]
    
    for guess in qudit_stabs:
         ops=[]
         print(' \n Currently checking if' ,guess,'is a stabiliser ')
         for j in range(n):
            temp_new=guess[j]
            
            if temp_new[0]=='X':
                x_number=temp_new[2]
                ops.append(np.around(xs[int(x_number)],6))
                 
                 
                 
                 
            elif temp_new[0]=='Z':
                z_number=temp_new[2]
                ops.append(np.around(zs[int(z_number)],6))
            elif temp_new[0]=='I':
                 ops.append(np.identity(d))
                 
         for k in range(n-1):
            ops[0]=np.kron(ops[0],ops[1])
            del(ops[1])
         #print(ops, '\n')
        # print('\n',ops[0],guess ,'\n')
         test=np.matmul(ops[0],input_state)
         if np.array_equal(np.around(test,2),np.around(input_state,2)):
             stabiliser_list.append(guess)
             stab_mats.append(ops[0])
             print(guess, 'is a stabiliser \n')
             
    
    #This just removes the multiple identities from the collection of stabilisers
    I_comp=np.identity(d**n,dtype=complex)
    qudit_stabiliser_matrices=[]
    checked_stabiliser_list=[]
    for i in range(len(stab_mats)):
        if np.array_equal(np.around(stab_mats[i]),I_comp):
            continue
        elif stabiliser_list[i].count('Z(0)')>0 or stabiliser_list[i].count('X(0)')>0 :
            continue
        else:
            qudit_stabiliser_matrices.append(stab_mats[i])
            checked_stabiliser_list.append(stabiliser_list[i])
    
    
    
    
    print(' \n \n \n \n \n \n The state:', np.around(input_state,3), 'is a qudit graph state with the stabilsers: \n ',checked_stabiliser_list,'\n\n')#, qudit_stabiliser_matrices)
   # for q in range(len(checked_stabiliser_list)):
    #    print(checked_stabiliser_list[q], '\n\n\n\n', qudit_stabiliser_matrices[q])
    #print(qudit_stabiliser_matrices)
    
    minimal_generator_set=[]
    for i in range(len(checked_stabiliser_list)):
        if checked_stabiliser_list[i][i]=='X(1)':
            minimal_generator_set.append(checked_stabiliser_list[i])
    print('The minimal generating set of the stabiliser group is:\n', minimal_generator_set)
    
    
    return checked_stabiliser_list, qudit_stabiliser_matrices, minimal_generator_set





#Function that creates an arbritrary qudit graph state vector 
#given the adjacency matrix, d and n
def create_qudit_graph_state_vector(d,n,A):
    init_state=[]
    init_string=str()
    w=cmath.exp(2*math.pi*1j/d)
     
    for t in range(d**n):
        init_state.append(1)
        
        if t<d:
            init_string=init_string+str(t)  
        
    input_strings=["".join(seq) for seq in itertools.product(init_string, repeat=n)]
    
    
    for i in range(n):
        for z in range(n):
            
            if z<=i:        #Only loop through upper triangle of adj matrix
                continue
            
            elif A[i][z]!=0:
                for f in range(len(input_strings)):
                    temp=input_strings[f]
                    j=int(temp[i])
                    k=int(temp[z])
                    w_temp=w**(j*k*A[i][z])
                    init_state[f]=init_state[f]*w_temp
            else:
                continue
    return np.around(init_state,6)


d=7
n=3

#State creation to test formalism###########################################################
A=np.array([[0,1,2],[2,0,1],[2,1,0]])
#A=np.array([[0,1,2,1],[2,0,1,0],[2,1,0,2],[1,0,2,0]])
#A=np.array([[0,1,2,1,0],[2,0,1,0,0],[2,1,0,2,0],[1,0,2,0,1], [0,0,0,1,0]])
#A=np.array([[0,2],[2,0]])
input_state=create_qudit_graph_state_vector(d,n,A)

# =============================================================================
# #input_state=[1+0j,1+0j,1+0j,-1+0j]
# #input_state=[1+0j,1+0j,1+0j,1+0j,1+0j,-1+0j,-1+0j,1+0j] #3 qubit line
# #                                                       #4 qubit star
# 
# p=-1/2+(math.sqrt(3)/2)*1j
# t=-1/2-(math.sqrt(3)/2)*1j
# 
# #input_state=[1+0j,1+0j,1+0j,1+0j,p,t,1+0j,t,p]  #2 qutrit 1 edge
# #input_state=[1+0j , 1+0j ,1+0j ,1+0j ,t ,p ,1+0j ,p ,t]   #2 qutrit 2 edges
# 
# =============================================================================
###################################################################################

stab_list, stab_matrices, min_gen_set=qudit_brute_force_stabiliser_check(input_state,d,n)
gp.get_qudit_graph_from_stabilisers(min_gen_set,d)




