#!/anaconda/bin/python
import os
import sys
import csv
import subprocess
from subprocess import Popen, PIPE

#os.chdir("/Users/dbauer/Research/DiagHam/FTI/src/Programs/FCI/")
os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))
filename = sys.argv[1]
p = int(sys.argv[2])

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

evals = []
if(filename == "-r"):
    p1 = subprocess.Popen(["ls","-t"],stdout=PIPE)
    p2 = subprocess.Popen(["head","-n 1"],stdout=PIPE,stdin=p1.stdout)
    #p3 = subprocess.Popen(["tail","-n 1"],stdout=PIPE,stdin=p2.stdout)
    filename = p2.communicate()[0].rstrip()
    print filename

f = open(filename,'r')
reader = csv.reader(f,delimiter=' ')
for row in reader:
    if(is_number(row[1])):
        evals.append(float(row[2]))
s = sorted(evals)

for i in range(int(p)+1):
    print s[i]

print str(s[p] - s[p-1])