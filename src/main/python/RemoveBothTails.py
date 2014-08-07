#!/usr/bin/python

import sys
import os

if len(sys.argv)>=4:
    filetype = sys.argv[1]
    sub_filename = sys.argv[2]
    output_filename = sys.argv[3]
else:
    print("usage: ./RemoveBothTails.py filetype sub_file output_sub_filename")
    print("or ")
    sys.exit(1)
################################################################################
def processls(d):
    result_d = {}
    i=0
    for range in d:
        ls = range.split('_')
        L = int(ls[1]) -int(ls[0])
        if i==0:
            ref_L = L
            ref_seq = d[range]
            ref_range = range
            i+=1 
            continue
        if L >= ref_L :
            result_d[range]=d[range]
        elif L < ref_L:
            result_d[ref_range]=ref_seq
            ref_L =L
            ref_seq = d[range]
            ref_range = range
    return result_d 
################################################################################

if (filetype == "fa" ):
    start_char = ">"
else:
    start_char = "@"    

sub_file = open(sub_filename,'r')
output = open(output_filename,'w')
intact_MS = open(output_filename+"_intact_MS",'w')
TF = 0
sub_dt ={}
while (True):
    line = sub_file.readline()
    if (line == ""):
        break
    if (line[0] == start_char):
        ls = (">" + line[1:]).strip().split('/')
        if ls[-1] == "ccs":
            output.write(">" + line[1:])
            TF = 1
        else:
            sub_name = '/'.join(ls[0:-1])+ '/'
            if not sub_dt.has_key(sub_name):
                sub_dt[sub_name] = {}
            sub_dt[sub_name][ls[-1]] = ""
            TF = 0
        continue
    if TF:
        output.write(line)
    else:
        sub_dt[sub_name][ls[-1]] = sub_dt[sub_name][ls[-1]] + line.strip()
    if (filetype == "fq"):
        line = sub_file.readline()  # skip quality lines
        if (line[0] != "+"):
            print "Err in LR fastq file format"
            exit(1)
        line = sub_file.readline()
            
sub_file.close()

for sub_name in sub_dt:
    if len(sub_dt[sub_name])>2:
        intact_MS.write( sub_name[1:] + '\t' + str(len(sub_dt[sub_name])) + '\n')
    if len(sub_dt[sub_name])==2:
        sub_dt[sub_name] = processls(sub_dt[sub_name])
    elif len(sub_dt[sub_name])>2:
        sub_dt[sub_name] = processls(sub_dt[sub_name])
        sub_dt[sub_name] = processls(sub_dt[sub_name])
    for range in sub_dt[sub_name]:
        output.write(sub_name+range+'\n')
        output.write(sub_dt[sub_name][range]+'\n')
output.close()
intact_MS.close()
