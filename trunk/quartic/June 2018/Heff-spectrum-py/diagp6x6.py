#!/opt/local/bin/python
import numpy as np
import scipy.linalg
import math
import sys

dim = int(sys.argv[1])

Hamiltonian = np.matrix(np.zeros(shape=(dim,dim)),dtype=complex)

f = open("2n=6_gs_evec" + str(dim) + ".dat",'w')
g = open("2n=6_evals" + str(dim) +".dat",'w')
h = open("2n=6_1es_evec.dat",'w')

for i in range(0,dim):
    for j in range(0,dim):
        if(i == j):
            Hamiltonian[i,i] = -(15.0/4.0) + (45.0/2.0)*(i+1) - (45.0/2.0)*(i+2)*(i+1) + 5.0*(i+3)*(i+2)*(i+1)
        elif(i == j - 4):
            Hamiltonian[i,j] =  -(15.0/4.0)*math.sqrt(j*(j-1)*(j-2)*(j-3)) + (3.0/2.0)*(j+1)*math.sqrt(j*(j-1)*(j-2))

for i in range(0,dim):
    for j in range(0,dim):
        if(i > j):
            Hamiltonian[i,j] = Hamiltonian[j,i]

eigenvalues,eigenvectors = scipy.linalg.eigh(Hamiltonian)

k = 0
for e in np.nditer(eigenvalues):
    g.write(str(e) + "\n")
    if(k < 10):
        print(e)
    k = k +1

for d in eigenvectors[0]:
    f.write(str(d.real) + "\n")
for d in eigenvectors[1]:
    h.write(str(d.real) + "\n")


h.close()
f.close()
g.close()