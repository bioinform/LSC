#!/usr/bin/python

import sys

if len(sys.argv)>=3:
    in_filename = sys.argv[1]
    out_filename = sys.argv[2]

else:
    print("Remove quality information and separate the reads from their names for performance reasons")
    print("usage: ./FASTQ2fa.py input_filename output_filename")
    print("or python FASTQ2fa.py input_filename output_filename")
    sys.exit(1)
################################################################################
infile = open(in_filename,'r')
outfile = open(out_filename,'w')
indexfile = open(out_filename+".readname",'w')

i=1
line_i = 0
for line in infile:
    if (line_i == 0):
        if line[0] !='@':
            print "Err: invalid LR fastq format"
            exit(1)
        outfile.write('>' + str(i) + '\n')
        indexfile.write(str(i) + '\t' + line[1:-1] + '\n')
        i+=1
        line_i = 1
    elif (line_i == 1):
        outfile.write(line)
        line_i = 2
    elif (line_i == 2):
        line_i = 3
    elif (line_i == 3):
        line_i = 0        


indexfile.close()
outfile.close()
infile.close()
################################################################################
    
