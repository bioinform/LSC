#!/usr/bin/python

import sys
import datetime

t0 = datetime.datetime.now()
################################################################################
def compress(seq):
    # Return values
    cps_seq = seq[0] # Compressed sequence
    pos_ls=[]
    len_ls = []
    n_count = 0
    
    # Loop carry values
    repeat_count=0
    ref_s = seq[0]
    i=0
    for s in seq:
        if not ref_s == s:
            if (ref_s == 'N'):
                n_count += 1
            if repeat_count>1:
                len_ls.append(str(repeat_count))
                pos_ls.append(str(i)) 
            cps_seq = cps_seq + s
            repeat_count=1
            ref_s = s
            i+=1
        else:
            repeat_count+=1
    if (ref_s == 'N'):
        n_count += 1
    if repeat_count>1:
        len_ls.append(str(repeat_count))
        pos_ls.append(str(i)) 
    return cps_seq,pos_ls,len_ls,n_count

################################################################################
MinNonN=40
MaxN=1
if len(sys.argv) >= 4:
    for opt in sys.argv[1:]:
        if opt[0]=='-':
            opt_ls = opt.split('=')
            if opt_ls[0]=="-MinNonN":
                MinNonN = int(opt_ls[1])
            elif opt_ls[0]=="-MaxN":
                MaxN = int(opt_ls[1])
    filetype = sys.argv[-3]
    inseq_filename =  sys.argv[-2]
    outseq_prefix = sys.argv[-1]
else:
    print("Perform homopoylmer compression on a fasta or fastq file and generate a cps (fasta) and idx file")
    print("usage: python compress.py [-MinNonN=39] [-MaxN=1] filetype inseq out_prefix")
    print("or ./compress.py [-MinNonN=39] [-MaxN=1] filetype inseq out_prefix")
    sys.exit(1)

################################################################################

################################################################################
if (filetype == "fa"):
    start_char = ">"
elif (filetype == "fq"):
    start_char = "@"    

inseq=open(inseq_filename,'r')
outseq=open(outseq_prefix+'cps','w')
idx = open(outseq_prefix+'idx','w')

# Process all the entries, one per iteration
while(True):
    # Read in the readname
    line = inseq.readline()
    if (line == ""):
        break
    readname = line[1:-1]
    if (line[0] != start_char):
        print line
        print "Err:  invalid file format: " + inseq_filename
        exit(1)
    
    # Read in the sequence
    line = inseq.readline()
    if (line == ""):
        break
    seq = line.strip().upper()
    
    # Compress, filter and output the sequence
    cps_seq, pos_ls, len_ls, n_count = compress(seq)
    if len(cps_seq)-n_count>=MinNonN and n_count<=MaxN:
        outseq.write(">"+readname+'\n' + cps_seq+'\n')
        idx.write(readname + "\t" + ','.join(pos_ls) + '\t' + ','.join(len_ls) + '\n')
    
    # Skip FASTQ quality lines
    if (filetype == "fq"):
        line = inseq.readline()
        line = inseq.readline()

inseq.close()
outseq.close()
idx.close()
print "finsish genome"

################################################################################

