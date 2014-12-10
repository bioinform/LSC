#!/usr/bin/python

import sys

if len(sys.argv) >= 2:
    filename = sys.argv[1]
else:
    print("usage: ./uniqseq2fasta.py SR_uniq.seq")
    print("or pythn SR_uniq.seq")
    sys.exit(1)
################################################################################
sr_uniq_file = open(filename ,'r')

i=1
for line in sr_uniq_file:
    ls = line.strip().split(" ")
    print ">" + str(i) + "_" + ls[-2]
    print ls[-1]
    i += 1
sr_uniq_file.close()
################################################################################
    
