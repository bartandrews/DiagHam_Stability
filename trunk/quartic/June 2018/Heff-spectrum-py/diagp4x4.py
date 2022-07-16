#!/opt/local/bin/python
import numpy as np
import scipy.linalg
import math
import sys

delta = 0
dim = int(sys.argv[1])

Hamiltonian = np.matrix(np.zeros(shape=(dim,dim)),dtype=complex)

f = open("2n=4_gs_evec" + str(dim) + ".dat",'w')
g = open("2n=4_evals" + str(dim) + ".dat",'w')
h = open("2n=4_1es_evec" +  str(dim) + ".dat",'w')
h2 = open("overlaps" +  str(dim) + ".dat",'w')


#quadratic hamiltonian
# for i in range(0,dim):
#     for j in range(0,dim):
#         if(i == j):
#             Hamiltonian[i,i] = Hamiltonian[i,i] + i + 1/2



#quartic hamiltonian
for i in range(0,dim):
    for j in range(0,dim):
        if(i == j):
            Hamiltonian[i,i] = (3/8.0) + 3.0/(4.0)*i + 3.0/(4.0)*i*i
        elif(i == j + 4):
            Hamiltonian[i,j] = (1/8.0)*math.sqrt((j+1)*(j+2)*(j+3)*(j+4))
        elif(i == j - 4):
            Hamiltonian[i,j] = (1/8.0)*math.sqrt(j*(j-1)*(j-2)*(j-3))

#quadratic hamiltonian
for i in range(0,dim):
    for j in range(0,dim):
        if(i == j):
            Hamiltonian[i,i] = Hamiltonian[i,i] + delta*(i + 1/2)



eigenvalues,eigenvectors = scipy.linalg.eigh(Hamiltonian)

k = 0
for e in np.nditer(eigenvalues):
    g.write(str(e) + "\n")
    if(k < 10):
        print(e)
    k = k +1
tracesum = 0
for (index,coef) in enumerate(eigenvectors[0]):
    tracesum = tracesum + np.conj(coef)*coef*index
print(tracesum)

h0 = 0
for (index,coef) in enumerate(eigenvectors[0]):
    h0 = h0 + (2*index+1)*np.conj(coef)*coef
print(h0)

for d in eigenvectors[0]:
    f.write(str(d.real) + "\n")
for d in eigenvectors[1]:
    h.write(str(d.real) + "\n")

for i in range(50):
    ev = eigenvectors[i]
    h2.write(str(abs(ev[i].real)) + "\n")

h2.close()
h.close()
f.close()
g.close()