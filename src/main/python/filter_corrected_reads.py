#!/usr/bin/python

import sys

if len(sys.argv)>=3:
    coverage_threshold = float(sys.argv[1]) # Threshold for percentage of the LR covered by SRs, stored in the corrected LR readname
    input_filename = sys.argv[2]
else:
    print("Post-runLSC script which can be used to filter the corrected long reads")
    print("usage: ./filter_corrected_reads.py coverage_threshold input_filename")
    print("or python filter_corrected_reads.py coverage_threshold input_filename")
    sys.exit(1)


input_file = open(input_filename,'r')

valid = 0
# For each read, filter it according to coverage 
for line in input_file:
    if ((line[0] == ">") or (line[0] == "@") ):
        fields = line.strip().split('|')
        coverage = float(fields[-1])
        if (coverage >=  coverage_threshold):
            print(line.strip())
            valid = 1
        else:
            valid = 0
    else:
        if (valid):
            print(line.strip())
            
input_file.close()
