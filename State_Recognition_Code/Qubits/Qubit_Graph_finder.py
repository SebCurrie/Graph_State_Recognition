import numpy as np
import itertools
import copy
import Qubit_Graph_methods as gp
import math
from scipy import sparse

#Functions that check an input state is a graph state by testing all combinations of stabilisers

#TO DO: CHANGE TO USE SCIPY SPARSE MATRICES

#Two qubit matrix definitions
hadamard=sparse.csc_matrix([[1,1],[1,-1]])
identity=sparse.csc_matrix([[1,0],[0,1]])
Z=sparse.csc_matrix([[1,0],[0,-1]])
X=sparse.csc_matrix([[0,1],[1,0]])




#Takes an input state and tests if its a graph state, if it is, returns the generators of the stabilser group
def brute_force_stab_check(input_coeffs,n):
    stabs=[]
    checked_stabiliser_list=[]
    stab_mats=[]
   
    #This generates a cartesian product of all possible combs of Z and I
    #doing like this allows us to just input n. We then copy the list
    #multiple times and insert X into a different position on each list 
    combo=itertools.product(['Z', 'I'], repeat=n-1)
    combos=[list(p) for p in combo]
 
    for i in range(n):
        stabs.extend(copy.deepcopy(combos))

    for i in range(n):
        for j in range(2**(n-1)):
           stabs[i*(2**(n-1))+j].insert(i,'X')
           
    print('The set of all possible stabilsers to be tested for n=',n,'is: \n', stabs, '\n')
    
    
    
    # Now given the list of all possible combinations we need to test, we generate a matrix for each
    #combination and multiply it by the input state vector. If its an eigenvector we append it to the list of 
    # checked stabilisers
    for guess in stabs:
        ops=[]
        print('Currently checking if' ,guess,'is a stabiliser \n')
        for j in range(n):
        
            if guess[j]=='X':
                ops.append(X)
            elif guess[j]=='Z':
                ops.append(Z)
            elif guess[j]=='I':
                ops.append(identity)
                
        for k in range(n-1):
            ops[0]=sparse.kron(ops[0],ops[1])
            del(ops[1])
            
        test=sparse.csc_matrix.dot(ops[0],input_coeffs)
        if np.array_equal(test,input_coeffs):
            checked_stabiliser_list.append(guess)
            stab_mats.append(ops[0])
    

    if len(stab_mats)!=n:
        print('This is not a graph state as there are not n generators. \n' )
        return _,_
    
    else:
        print('This is a graph state with stabilsers generators:', checked_stabiliser_list,'\n')#,stab_mats)
        return checked_stabiliser_list,stab_mats
    
    
    
    
    
    
#FUNCTION NOT WORKING YET
#This function doesnt work properly yet, returns >2**n stabilisers
# =============================================================================
# def all_stabs_from_generators(generator_list, generator_matrices):
#     #Now we need to double check that every possible combination of these are also stabilsers(they should be).  
#     stab_mats=generator_matrices
#     first=True       
#     for i in range(n):
#         for j in range(n):
#             stab_comb=np.matmul(stab_mats[i],stab_mats[j])
#     
#             if i==j and first==True and  np.array_equal(np.matmul(stab_comb,input_coeffs),input_coeffs): 
#                 first=False                                                                                 #Avoids ending up with n identites in the list
#                 stab_mats.append(stab_comb) 
#                 #checked_stabiliser_list[i].append(checked_stabiliser_list[j])
#                 
#             elif i!=j and np.array_equal(np.matmul(stab_comb,input_coeffs),input_coeffs):#
#                 stab_mats.append(stab_comb) 
#                # checked_stabiliser_list[i].append(checked_stabiliser_list[j])
#                 
#             elif i!=j:
#                 print('The combination of stabilisers is not a stabilser??')
#     
#     
#     
#     #Removing any other duplicates
#     unique_stabs = []
#     for arr in stab_mats:
#         if not any(np.array_equal(arr, unique_arr) for unique_arr in unique_stabs):
#             unique_stabs.append(arr)
#     
#     
#     for i in range(len(unique_stabs)):
#         print(unique_stabs[i],'\n',checked_stabiliser_list, '\n')
#     
#     if len(unique_stabs)!= 2**n:
#         print('Not a graph state there should be 2**n=',2**n,'stabilsers but found', len(stab_mats), 'stabilisers')
#         
# 
# 
#     return all_stabilisers
# 
# =============================================================================




#HARD CODING INPUT TO TEST FUNCTIONS===============================================================

n=3
#input_coeffs=[1,1,1,-1]
input_coeffs=[1,1,1,-1,1,1,-1,1]

#================================================================================================


generator_list, generator_matrices=brute_force_stab_check(input_coeffs,n)
gp.get_graph_from_stabilisers(generator_list)



