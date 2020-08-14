import itertools
import Qudit_graph_methods as gp
import math
import numpy as np
import copy
import cmath
from scipy import sparse 

#WARNING: THIS CODE IS EXTREMELY COMPUTATIONALLY EXPENSIVE ABOVE 5 QUBITS, AND A FAR MORE EFFICIENT ALGORITHM IS IMPLEMENTED IN QUDIT HYPERGRAPH FINDER
#THIS CODE SHOULD ONLY REALLY BE USED AS A CHECK FOR THE OUTPUT OF THAT.


#This file contains the functions neccessary to perform a brute force stabiliser check on a given input state
#This will check if the state is a graph state of qudits, but does not have the functionality to 
#recognise hypergraphs of qudits.



#Generate the "sub" matrices before combining to create elements of the pauli group
#Will return d matrices for b=0,1,2...d-1
def generate_generalised_pauli_matrices(d):
    w=cmath.exp(2*math.pi*1j/d)
    
    zs=[]
    for l in range(d):
        z=np.array(np.zeros((d)),dtype=complex)
        
        for t in range(d):
            z[t]=w**(l*t)
        z_sparse=sparse.diags(z,offsets=0)
        zs.append(z_sparse)
        
    
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
        x_sparse=sparse.csc_matrix(x)
        print(x_sparse)
        xs.append(x_sparse)
    print(xs)
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
    #Now we are in the qudit case there are multiple Z matrices to check and we must check all possible combs
#   i.e X(0)Z(0), X(0)Z(1), X(1),Z(0), X(1)Z(1)
#       Z(0)X(0), Z(1)X(0),  Z(0)X(1)  , Z(1)X(1)
    
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
            ops[0]=sparse.kron(ops[0],ops[1])
            del(ops[1])
         #print(ops, '\n')
        # print('\n',ops[0],guess ,'\n')
      
         test=sparse.csc_matrix.dot(ops[0],input_state)
        # print(type(test))
         if np.array_equal(np.around(test,2),np.around(input_state,2)):
             stabiliser_list.append(guess)
             stab_mats.append(ops[0])
             print(guess, 'is a stabiliser \n')
             
    
    #This just removes the multiple identities from the collection of stabilisers
    qudit_stabiliser_matrices=[]
    checked_stabiliser_list=[]
    for i in range(len(stab_mats)):
        if stabiliser_list[i].count('Z(0)')>0 or stabiliser_list[i].count('X(0)')>0: #or stabiliser_list[i].count('I')>(n-1):
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
    
    gp.nx_get_qudit_graph_from_stabilisers(minimal_generator_set,d)
    
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





#SECTION FOR USER INPUT TO CALL AND TEST FUNCTIONS=================================================
    
d=7
n=5

#A=np.array([[0,1,2],[2,0,1],[2,1,0]])
#A=np.array([[0,1,2,1],[2,0,1,0],[2,1,0,2],[1,0,2,0]])
A=np.array([[0,1,2,1,0],[1,0,1,0,0],[2,1,0,2,0],[1,0,2,0,1], [0,0,0,1,0]])
#A=np.array([[0,1,2,1,0,0],[2,0,1,0,0,0],[2,1,0,2,0,0],[1,0,2,0,1,0], [0,0,0,1,0,1],[0,0,0,1,1,0]])
#A=np.array([[0,2],[2,0]])



input_state=create_qudit_graph_state_vector(d,n,A)
stab_list, stab_matrices, min_gen_set=qudit_brute_force_stabiliser_check(input_state,d,n)




