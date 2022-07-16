#! /usr/bin/env python

# -*- coding: utf-8 -*-

import os,glob
import sys


file1=sys.argv[1]

print file1
f1=open(file1,"r")


bib_lines=f1.readlines()

bib_out_lines=[]
for line1 in bib_lines:
    if line1.startswith("file"):
        print line1
        p1=line1.find("Users")
        file_locn="/"+line1[p1:len(line1)-7]
        file_name=os.path.split(file_locn)[1]
        print file_name
        if not(os.path.isdir("./BibDesk-Papers")):
        	print "dir doesnt exist"
        	os.mkdir("./BibDesk-Papers")
        	#break
        cmd1='cp "'+ file_locn + '" ./BibDesk-Papers/'
        print cmd1
        output2=os.popen(cmd1).read()
        print output2
        #break
        p2=line1.find("localhost")
        dir_name=os.getcwd()
        line2="Local-Url = {file://localhost" + dir_name.replace(" ","%20",)+ "/BibDesk-Papers/"+file_name.replace(" ","%20")+ "},\n"
        print line2
        #bib_out_lines.append(line2)
        line3="file = {:."+   "/BibDesk-Papers/"+file_name+ ":pdf},\n"
        print line3
        bib_out_lines.append(line3)
        #break
    else:
        bib_out_lines.append(line1)


f2=open("out.bib","w")
f2.writelines(bib_out_lines)
f2.close()
#print bib_out_lines
#continue



