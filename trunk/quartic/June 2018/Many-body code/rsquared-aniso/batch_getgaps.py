#!/opt/local/bin/python

import sys
import os
import csv
import os.path

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))
batchfile_name = sys.argv[1]
outfile_name = sys.argv[2]
batchfile = open(batchfile_name,'r')

reader = csv.reader(batchfile,delimiter=',')
i = 1
data = []
for row in reader:
    evals = []
    np = int(row[0])
    flux = 1
    lx = int(row[1])
    ly = int(row[2])
    X = int(row[3])
    Y = int(row[4])
    t2 = float(row[5])
    gx = 0.0
    gy = float(row[6])
    if(t2==0):
        t2=int(0)
    t2 = '{:.2f}'.format(t2)
    t2 = float(t2)
    if(t2==0):
        t2=int(0)
    if(gy==0):
        gy=int(0)
    if(flux > 1):
        multiband_str = "_multiband_"
    else:
        multiband_str = "_"
    filename = ("fermions_hofstadter_X_" + str(X) + "_Y_" + str(Y) + "_q_" + str(flux) + "_n_" + str(np) + "_x"
                "_" + str(lx) + "_y_" + str(ly) + "_t2_" + str(t2) + "_t3_0_alpha_" + str(gy) +"_gx_0_gy_0.dat")
    #print(filename)
    if(os.path.isfile(filename)):
        print("Opening " + filename)
        myfile = open(filename,'r')
        gap_reader = csv.reader(myfile,delimiter=' ')
        for row in gap_reader:
            if(is_number(row[1])):
                evals.append(float(row[2]))
        s = sorted(evals)
        print(s[4])
        print(s[3])
        print(s[2])
        delta = float(s[3]-s[2])
        lx = float(lx)
        ly = float(ly)
        X = float(X)
        Y = float(Y)
        print(((X*Y)**2)*delta)
        totX = lx*X
        totY = ly*Y
        aniso = totY/totX
        data.append([int(np),int(X*Y),float(t2),str(aniso),str(gy),((X*Y)**2)*float(delta)])
        myfile.close()
    else:
        print("File " + filename + " not found.")
batchfile.close()

with open(outfile_name,'w') as outfile:
    cwrite = csv.writer(outfile,delimiter=',')
    for d in data:
        cwrite.writerow(d)

         
