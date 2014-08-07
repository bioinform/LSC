#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 2:
    filename = sys.argv[1]
else:
    print("usage: ./uniqseq2fasta.py SR_uniq.seq")
    print("or pythn SR_uniq.seq")
    sys.exit(1)
################################################################################
file = open(filename ,'r')

i=1
for line in file:
    ls = line.strip().split(" ")
    print ">" + str(i) + "_" + ls[-2]
    print ls[-1]
    i += 1
file.close()
################################################################################
    
