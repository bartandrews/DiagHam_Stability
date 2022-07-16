#!/home/bart/miniconda3/bin/python

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
batchfile = open(batchfile_name, 'r')

reader = csv.reader(batchfile, delimiter=',')
i = 1
data = []

data.append(["Np", "UCarea", "t2", "aniso", "twist", "gap"])

for row in reader:
    evals = []
    np, lx, ly, X, Y, t2, gy = int(row[0]), int(row[1]), int(row[2]), int(row[3]), int(row[4]), float(row[5]), float(row[6])
    flux = 1
    if t2 == 0:
        t2 = int(0)
    t2 = '{:.2f}'.format(t2)
    t2 = float(t2)
    if t2 == 0:
        t2 = int(0)
    if gy == 0:
        gy = int(0)
    if flux > 1:
        multiband_str = "_multiband_"
    else:
        multiband_str = "_"
    filename = ("bosons_hofstadter_X_" + str(X) + "_Y_" + str(Y) + "_q_" + str(flux) + "_n_" + str(np) + "_x"
                "_" + str(lx) + "_y_" + str(ly) + "_t2_" + str(t2) + "_t3_0_alpha_1_u_1_gx_0_gy_" +
                str(gy) + ".dat")

    if os.path.isfile(filename):
        print("Opening " + filename)
        myfile = open(filename, 'r')
        gap_reader = csv.reader(myfile, delimiter=' ')
        for row in gap_reader:
            if is_number(row[1]):
                evals.append(float(row[2]))
        s = sorted(evals)
        delta = float(s[1]-s[0])
        lx = float(lx)
        ly = float(ly)
        X = float(X)
        Y = float(Y)
        print("Gap:")
        print(((X*Y))*delta)
        totX = lx*X
        totY = ly*Y
        aniso = totY/totX
        data.append([int(np), int(X*Y), float(t2), str(aniso), str(gy), float(X*Y*delta)])
        myfile.close()
    else:
        print("File " + filename + " not found.")
batchfile.close()

with open(outfile_name, 'w') as outfile:
    cwrite = csv.writer(outfile, delimiter=',')
    for d in data:
        cwrite.writerow(d)

         
