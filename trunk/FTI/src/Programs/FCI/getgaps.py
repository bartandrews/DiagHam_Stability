import os
import os.path
import sys
import csv

#os.chdir("/Users/dbauer/Research/DiagHam/FTI/src/Programs/FCI/")
os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))
outfile = open("gapsVSN_tau18.csv",'w')
outfile_f = open("gapsVSNf.csv",'w')
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
evals = []
for i in range(2,14):
    evals = []
    size = i
    filename = "bosons_hofstadter_X_" + str(i) + "_Y_" + str(i) + "_q_1_sym_n_8_x_4_y_4_gx_0_gy_0.dat"
    myfile = open(filename,'r')
    reader = csv.reader(myfile,delimiter=' ')
    for row in reader:
        if(is_number(row[1])):
            evals.append(float(row[2]))
    s = sorted(evals)
    k = 1
    print i
    #for p in s:
    #    if k < 4:
    #        print p
    #    k = k+1
    print str((s[2]-s[1]))
    outfile.write(str(i*i)+","+str(s[2]-s[1])+"\n")

for j in range (1,8):
    evals = []
    filename = "fermions_hofstadter_X_" + str(j*6) + "_Y_" + str(j) + "_q_1_n_8_x_2_y_12_gx_0_gy_0.dat"
    myfile = open(filename,'r')
    reader = csv.reader(myfile,delimiter=' ')
    for row in reader:
        if(is_number(row[1])):
            evals.append(float(row[2]))
    s = sorted(evals)
    #for q in s:
        #print q
    print j
    print str((s[3]-s[2]))
    outfile_f.write(str(6*j*j)+","+str(s[3]-s[2])+"\n")
