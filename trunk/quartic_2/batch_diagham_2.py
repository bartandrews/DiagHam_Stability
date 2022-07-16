#!/home/bart/miniconda3/bin/python

import sys
import os
import csv
import subprocess
import time


os.chdir(os.path.dirname(os.path.abspath(sys.argv[1])))
print(os.getcwd())
batchfile_name = sys.argv[1]
FNULL = open(os.devnull, 'w')
batchfile = open(batchfile_name,'r')

reader = csv.reader(batchfile,delimiter=',')
i = 1
for row in reader:
    cmd_string = "/home/bart/DiagHam_Stability/trunk/build/FTI/src/Programs/FCI/FCIHofstadterModel --boson --lanczos-precision 1e-10 -n 2 -m 30000 --use-lapack "#"--three-body-potential 1.0 "
    # cmd_string.append("-p " + str(row[0]))
    # cmd_string.append("-x " + str(row[1]))
    # cmd_string.append("-y " + str(row[2]))
    # cmd_string.append("-X " + str(row[3]))
    # cmd_string.append("-Y " + str(row[4]))
    # cmd_string.append("-t2 " + str(row[4]))
    cmd_string = cmd_string + "-p " + str(row[0]) + " "
    cmd_string = cmd_string + "-x " + str(row[1]) + " "
    cmd_string = cmd_string + "-y " + str(row[2]) + " "
    cmd_string = cmd_string + "-X " + str(row[3]) + " "
    cmd_string = cmd_string + "-Y " + str(row[4]) + " "
    cmd_string = cmd_string + "--t2 " + str(row[5]) + " "
    cmd_string = cmd_string + "--gamma-y " + str(row[6]) + " "
#    cmd_string = cmd_string + "-t3 " + str(row[6]) + " "
    #cmd_string = cmd_string + "--eigenstate" + str(int(row[4])*int(row[5]))
    # if(int(row[1]) >= 2):
    #     cmd_string = cmd_string + "--multi-band"

    print("Running " + str(cmd_string))
    t1 = time.time()
    # os.system(cmd_string)  # execute the command
    process = subprocess.Popen(cmd_string,stdout=FNULL, stderr=FNULL,shell=True)
    process.wait()
    t2 = time.time()
    print("Run " + str(i) + " complete. Time " + str(t2-t1))
    i = i + 1


