import numpy as np
import itertools
import copy
import Graph_Printers as gp
import math

hadamard=np.array([[1,1],[1,-1]])
I=np.array([[1,0],[0,1]])
Z=np.array([[1,0],[0,-1]])
X=np.array([[0,1],[1,0]])


input_coeffs=[1,1,1,-1,-1,-1,-1,-1]

a=np.kron(X,I)
a=np.kron(a,I)

b=np.kron(I,Z)
b=np.kron(b,I)

CZ=np.array([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,-1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,-1]])


#out=np.matmul(a,b)
out=np.matmul(a,CZ)

state=np.matmul(out,input_coeffs)

print(out,state)