#!/opt/local/bin/python

import sys
import os
import csv
import subprocess
import time


os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))
batchfile_name = sys.argv[1]
FNULL = open(os.devnull, 'w')
batchfile = open(batchfile_name,'r')

reader = csv.reader(batchfile,delimiter=',')
i = 1
for row in reader:
    cmd_string = "./FCIHofstadterModel --boson --flat-band --lanczos-precision 10e-9 -n 2 -m 4000 --use-lapack "#"--three-body-potential 1.0 "
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
    process = subprocess.Popen(cmd_string,stdout=FNULL, stderr=FNULL,shell=True)
    process.wait()
    t2 = time.time()
    times=[]
    times.append(t2-t1)
    #avg_t = sum(times)/len(times)
    print("Run " + str(i) + " complete. Time " + str(times[-1]))
    i = i + 1


