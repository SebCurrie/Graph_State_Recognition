# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 16:58:42 2020

@author: mr19164
"""


 qudito_stabs=[]
    for i in range(len(qudit_stabs)):
        temp=qudit_stabs[i]
        no_zs=temp.count('Z(b)')
        print('no z', no_zs)
        count=0
        for k in range(n):
            if temp[k]=='Z(b)':
                count+=1
                
                for j in range(d):
                    zs=[]
                    qudit_stab='Z('+str(j)+')'
                    zs.append(qudit_stab)
                
                
                for i in range(d):
                    temp[k]=zs[i]
                    qudit_stabs.append(copy.deepcopy(temp))
                    if count==no_zs:
                        qudito_stabs.append(copy.deepcopy(temp))
        
            elif temp[k]=='I':
                qudito_stabs.append(copy.deepcopy(temp))
                continue
            else:
                continue
                     

# =============================================================================
# 
#     #Now we have our elements X(a),Z(b), W(c) for all a,b,c in Fp
#     #Now all combinations of these form the generators of the generalised pauli group
#     # At the moment i dont think this is the MINIMAL generating set
#     gen_pauli_generators={}
#     for i in range(d):
#         for j in range(d):
#             for k in range(d):
#                 
#                 gen_pauli_generator=[]
#                 gen_pauli_generator=ws[i]*np.matmul(xs[j],zs[k])
#                 gen_pauli_generators[str(i)+str(j)+str(k)]=gen_pauli_generator
# 
#                 #print(gen_pauli_generator,'\n',i,j,k)
#     print('The generalised pauli group generators are:')            
#     for key in gen_pauli_generators:
#         print(gen_pauli_generators[key],'\n')
#     return gen_pauli_generators
# 
# =============================================================================


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
                for j in range(d):
                    qudit_stab='X('+str(j)+')'
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
        # print(' \n Currently checking if' ,guess,'is a stabiliser ')
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
    
    
    
    
    print(' \n \n \n \n \n \n The state:', input_state, 'is a qudit graph state with the stabilsers: \n ',checked_stabiliser_list,'\n\n')#, qudit_stabiliser_matrices)
   # for q in range(len(checked_stabiliser_list)):
    #    print(checked_stabiliser_list[q], '\n\n\n\n', qudit_stabiliser_matrices[q])
    #print(qudit_stabiliser_matrices)
    
    minimal_generator_set=[]
    for i in range(len(checked_stabiliser_list)):
        if checked_stabiliser_list[i][math.floor(i/(d-1))]=='X(1)':
            minimal_generator_set.append(checked_stabiliser_list[i])
    print('The minimal generating set of the stabiliser group is:', minimal_generator_set)
    
    
    return checked_stabiliser_list, qudit_stabiliser_matrices, minimal_generator_set

