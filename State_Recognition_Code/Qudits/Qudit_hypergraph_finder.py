import numpy as np
import itertools
import Qudit_hypergraph_methods as gp
import math
import cmath

#Two qubit matrix definitions
hadamard=np.array([[1,1],[1,-1]])
identity=np.array([[1,0],[0,1]])
Z=np.array([[1,0],[0,-1]])
X=np.array([[0,1],[1,0]])



#Need to write a function that tests whether RUN state under QFT
# =============================================================================
# def is_state_LEQ_REW(state):
#     no_zero_coeffs=state.count(0)
#     checks=itertools.product([hadamard, identity], repeat=n)
#     checks=[list(p) for p in checks]
#     
#     for i in range(2**n):
#            temp=checks[i]
#            if no_zero_coeffs==0:
#                REW_TEST(state)
#                print('input state:',input_coeffs,'is already a REW state')
#                return input_coeffs,_
#            else:
#                for k in range(n-1):
#                    temp[0]=np.kron(temp[0],temp[1])
#                    del(temp[1])
#                LU_State=np.matmul(temp[0],input_coeffs)
#         
#                if REW_TEST(LU_State):
#                     print('State',input_coeffs,' is LU to the REW state',LU_State,' under the action of \n',temp[0],'\n')
#                     return LU_State, temp[0]
#     
# =============================================================================





#First we need to test whether the state we have is in a "Roots of Unity" (RUN) state.
#This can be thought of as a generalisation of REW states.
#If a state is not a RUN state it cannot be a qudit hypergraph state
#Function takes a state vector and the dimension of the qudits

def IS_RUN_State(state,d):
    
    roots=[]
    w=cmath.exp(2*math.pi*1j/d)
    for i in range(d):
        roots.append(np.around(w**i,4))  
    for i in range (len(state)):
        if np.around(state[i],4) not in roots and np.around(state[i],4)*(-1) not in roots:
            RUN=False
            print('Not a RUN state!!')
            break
        else:
            RUN=True
            continue
    
    if RUN:
        print('This state is a Roots of Unity (RUN) State.')
    
    return RUN
    






#Takes the coefficients and strings of a RUN state, 
#and returns the operations to take you to the associated hypergraph
#The arguements are the coefficients of the state vector, the strings of the state vector, dimension and number  of qudits
#This implements a generalisation of the algorithm in Rossi 2013
#Full detail can be found in the README
def qudit_hypergraph_from_RUN(coeffs,strings,d,n):
    ops=[]
    w=cmath.exp(2*math.pi*1j/d)
    string_no_zero=str()
    
    #used later to get the indices of non zero entries in a string
    for i in range(d-1):
        string_no_zero=string_no_zero+str(i+1)
    

    for t in range(n):  #First we are going to check for terms with one non-zero, then 2 ....
        for k in range(len(coeffs)):                        #Loop through every element of our state vector
            
            if np.around(coeffs[k],4)!=np.around(w,4) and strings[k].count('1')==t+1:      #If coefficient is not the first root of unity, and the string has t 1's in it
            
                weight=0
                test_coeff=coeffs[k] #Saves us having to skip this term in the loop later
                
                #The weight of the edge is the amount of times you have to apply CZ^-1 before getting back to w
                #This loops until the state coeff = w, and counts the no of appilcations of CZ^-1 required
                while np.around(test_coeff,7)!=np.around(w,7):                     
                    weight=weight+1
                    if weight>d:
                        break
                    
                    for g in range(t+1):
                        test_coeff=test_coeff/w**(int(strings[k][g])*weight)
                        
                
                #Now we check which qubits have the non-zero 'excitations'.
                indices = [i for i, x in enumerate(strings[k]) if x in string_no_zero]     #get the position of those 1's within the string. ie -1100 has ones at indexes 0,1
                ops.append(indices)
                
                
                #Now we need to perform (CZ^-1)**weight to every element with non-zeros at those indices
                for z in range(len(coeffs)):   
                    
                    #finds where the state vector element being checked has non zeros
                    to_be_tested_indices=[i for i, x in enumerate(strings[z]) if x in string_no_zero]      
           
                    #if the element has excitations in the same place we must apply the approptiate CZ^-1
                    if set(indices).issubset(set(to_be_tested_indices)):        
                        
                        #Calculates the appropriate power of w to divide by in order to apply CZ^-1
                        exponent=weight
                        for p in indices:
                            exponent=exponent*int(strings[z][p])
                        coeffs[z]=coeffs[z]/w**exponent
                    
                    else:
                        continue
                        
    
    
    #The space of RUN states is larger than the space of qudit hypergraphs (Xiong 2018)
    # There are states that are RUN states that will not be qudit hypergraphs, and if after 
    # the application of this algorithm, the state vector != Gphase* |+>^n then it is not a graph state 
    for i in coeffs:
        if np.around(i,4)!=np.around(w,4):
            raise Exception('The input RUN state did not correspond to a qudit hypergraph')
            
    
    if len(ops)!=0:            
        print('The edges in the graph are:',ops,'\n')
        gp.get_hypergraph_from_edges(ops,n)
    else:
        print('This is the empty hypergraph, will not display properly.')
  
    return ops
    
    

def generate_strings(n,d,input_coeffs):
    # #This is all possible combinations of bit strings
    # #useful later where we need index positions of states with one 1, two 1's etc
    init_string=str()
    for i in range(d):
        init_string=init_string+str(i)

    input_strings=["".join(seq) for seq in itertools.product(init_string, repeat=n)]
    if len(input_coeffs)!=len(input_strings):
        print('WARNING: Coefficient array and String array have different lengths')
        
    return input_strings



#User input to specify the state, dimensions and number of qubits
n=3
d=3
w=cmath.exp(2*math.pi*1j/d)
input_coeffs=[w,w,w,w,w,w,w,w,w,w**2,w**2,w**2,w**3,w**5,w**7,w**4,w**8,w**12,w**3,w**3,w**3,w**5,w**9,w**13,w**7,w**15,w**23]


#Function calls
state_strings=generate_strings(n,d,input_coeffs)
IS_RUN_State(input_coeffs,3)
qudit_hypergraph_from_RUN(input_coeffs,state_strings,d,n)


