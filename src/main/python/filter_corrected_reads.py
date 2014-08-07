#!/usr/bin/python

import sys
import os


coverage_threshold  =  float(sys.argv[1])
input_filename = sys.argv[2]

input_file = open(input_filename,'r')

valid = 0
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
