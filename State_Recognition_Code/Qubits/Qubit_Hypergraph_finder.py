import numpy as np
import itertools
import copy
import Qubit_Hypergraph_methods as gp
import math



#This is a python implementation of the algorithm in doi:10.1088/1367-2630/15/11/113022
#which identifies a qubit hypergraph state directly from the state vector
#meaning you dont have to do an expensive stabiliser check.

#A generalisation of this alogrithm to hypergraphs of qudits is available.



#Two qubit matrix definitions
hadamard=np.array([[1,1],[1,-1]])
identity=np.array([[1,0],[0,1]])
Z=np.array([[1,0],[0,-1]])
X=np.array([[0,1],[1,0]])




#Need to write a function that tests whether state is equivalent to a REW state under hadamards
#This function needs to test all possible combinations of Identity tensor hadamard
# III HII IHI ....HHH etc

def is_state_LEQ_REW(state):
    no_zero_coeffs=state.count(0)
    checks=itertools.product([hadamard, identity], repeat=n)
    checks=[list(p) for p in checks]
    
    for i in range(2**n):
           temp=checks[i]
           if no_zero_coeffs==0:
               REW_TEST(state)
               print('input state:',input_coeffs,'is already a REW state')
               return input_coeffs,_
           else:
               for k in range(n-1):
                   temp[0]=np.kron(temp[0],temp[1])
                   del(temp[1])
               LU_State=np.matmul(temp[0],input_coeffs)
        
               if REW_TEST(LU_State):
                    print('State',input_coeffs,' is LU to the REW state',LU_State,' under the action of \n',temp[0],'\n')
                    return LU_State, temp[0]
            
        
        
#Test whether a state is indeed a Real Equally Weighted state
def REW_TEST(state):
 
    for i in range (len(state)):
        if abs(state[0])!=abs(state[i]):
            REW=False
            break
        else:
            REW=True
            continue
        
    return REW
    


#Takes the coefficients and strings of a REW state, 
#and returns the operations to take you to the associated hypergraph         
def Graph_From_REW(coeffs,strings,n):
    ops=[]

    for t in range(n):  #First we are going to check for terms with one 1, then 2 ....
        
        for k in range(len(coeffs)):                        #Loop through every element of our bit string array
            if coeffs[k]<0 and strings[k].count('1')==t+1:      #If coefficient is negative, and the string has t 1's in it
                
                indices = [i for i, x in enumerate(strings[k]) if x == "1"]     #get the position of those 1's within the string. ie -1100 has ones at indexes 0,1
                ops.append(indices)
                
                for z in range(len(coeffs)):                        #Now we need to apply a negative to every element with this many 1's at those indices
                    to_be_tested_indices=[i for i, x in enumerate(strings[z]) if x == "1"]      #finds where the element has 1's
                    
                    if set(indices).issubset(set(to_be_tested_indices)):        #if the element has 1's in the same place we must apply a negative
                        coeffs[z]=-coeffs[z]
     
    if len(ops)!=0:            
        gp.get_hypergraph_from_edges(ops,n)
    else:
        print('This is the empty hypergraph, will not display properly.')
    
    gp.get_hypergraph_stabilisers(ops,n)
    
    return ops
    


#This generates all possible combinations of bit strings
#useful later where we need index positions of states with one 1, two 1's etc 
def generate_strings(n,input_coeffs):

    input_strings=["".join(seq) for seq in itertools.product("01", repeat=n)]
    

    if len(input_coeffs)!=len(input_strings):
        print('WARNING: Coefficient array and String array have different lengths')    


    return input_strings




#Section for hard coding input to test functions==========================================
    
#Non-REW input
#input_coeffs=[1,0,0,0,0,0,1,0] 

#REW input
n=3
input_coeffs=[1,1,1,-1,-1,-1,-1,-1]
#input_coeffs=[1,1,1,-1,1,1,-1,1]
#input_coeffs=[1,1,1,1,1,1,1,1]

#==========================================================================================

input_strings=generate_strings(n,input_coeffs)
state_coeffs,LU=is_state_LEQ_REW(input_coeffs)
edges=Graph_From_REW(state_coeffs, input_strings,n)
