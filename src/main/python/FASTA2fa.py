#!/usr/bin/python

import sys
import os

if len(sys.argv)>=3:
    in_filename = sys.argv[1]
    out_filename = sys.argv[2]

else:
    print("usage: ./FASTA2fa.py input_filename output_filename")
    print("or ptyhon FASTA2fa.py input_filename output_filename")
    sys.exit(1)
################################################################################
infile = open(in_filename,'r')
outfile = open(out_filename,'w')
indexfile = open(out_filename+".readname",'w')

i=1
for line in infile:
    if line[0]=='>':
         if i>1:
             outfile.write('\n')
         outfile.write('>'+str(i)+'\n')
         indexfile.write(str(i)+'\t'+line[1:-1]+'\n')
         i+=1
    else:
         outfile.write(line.strip())
outfile.write('\n')

indexfile.close()
outfile.close()
infile.close()
################################################################################
    
