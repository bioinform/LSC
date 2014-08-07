#!/usr/bin/python

import sys
import os
import re
import numpy

if len(sys.argv) >= 5:
    LR_filename = sys.argv[1]
    SR_filename = sys.argv[2]
    sam_filename = sys.argv[3]     #####Important Note: Input sam file is expected to be SR sorted
    nav_filename = sys.argv[4]
    one_line_per_alignment = (sys.argv[4][0] == "T")  # T: TRUE, F: FALSE
else:
    print("usage: python samParser.py LR.fa.cps SR.fa.cps sam_file nav_output_filename")
    print("or ./samParser.py LR.fa.cps SR.fa.cps sam_file nav_output_filename")
    sys.exit(1)

################################################################################

rev_cmplmnt_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 
                   'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
rev_cmplmnt_bases = rev_cmplmnt_map.keys()
def reverse_complement(seq):
    
    output_seq = ''
    len_seq = len(seq)
    for i in range(len_seq):
        if (seq[len_seq - 1 - i] in rev_cmplmnt_bases):
            output_seq += rev_cmplmnt_map[seq[len_seq - 1 - i]]
        else:
            print "Err: Unexpected base in short read sequence: " + seq
            output_seq += seq[len_seq - 1 - i]
        
    return output_seq
    

def get_SR_sequence(SR_file, SR_idx_file, SR_seq_readname):
    read_name = "invalid_read_name"
    while (read_name != SR_seq_readname):
        read_name = (SR_file.readline().strip())
        if (not read_name):
            print "Err: HC SR sequence not found."
            exit(1) 
        if (read_name[0] != '>'):
            continue    # unexpected string
        read_name = read_name[1:]
        SR_seq = SR_file.readline().strip().upper()
        SR_idx_seq = SR_idx_file.readline().strip().split('\t')
        if(SR_idx_seq[0] != read_name):
            print "Err: SR.fa.cps and SR.fa.idx do not match."
            exit(1) 
        SR_idx_seq = '\t'.join(SR_idx_seq[1:])
    return [SR_seq, SR_idx_seq]
    
LR_file = open(LR_filename,'r')
LR_seq = {}
line_num = 1
for line in LR_file:
    
    if (line_num):
        if (line[0] != '>'):
             continue    # unexpected string
        read_name = (line.strip())[1:]
    else:
        LR_seq[read_name] = line.strip().upper()
    line_num = 1 - line_num
LR_file.close()

SR_file = open(SR_filename + ".cps",'r')
SR_idx_file = open(SR_filename + ".idx",'r')
sam_file = open(sam_filename, 'r')
nav_file = open(nav_filename, 'w')

SR_seq_readname = "invalid_read_name"
SR_seq = ""
SR_seq_rvs_cmplmnt = ""
for line in sam_file:
    
    if (line[0] == '@'):
        continue
    
    line_fields = line.strip().split('\t')
    cigar = line_fields[5]
    if ((cigar == '*') or (cigar == '.')):
        continue
    
    SR_name = line_fields[0]
    if (SR_name != SR_seq_readname):
        [SR_seq, SR_idx_seq] = get_SR_sequence(SR_file, SR_idx_file, SR_name)
        SR_seq_rvs_cmplmnt = reverse_complement(SR_seq)
        SR_seq_readname = SR_name
    
    if (int(line_fields[1]) & 0x10):     # Check if seq is reversed complement
        line_fields[3] = '-' + line_fields[3]
    else:
        line_fields[3] = '+' + line_fields[3]
    
    align_list = [','.join([line_fields[2], line_fields[3], line_fields[5], str(0)])]
    
    if (not one_line_per_alignment):   # BWA reports all alignment per read in one line
        multi_align_str = ','.join([line_fields[2], line_fields[3], line_fields[5], str(0)]) + ';'
        for fields_idx in range(11, len(line_fields)):
            if (line_fields[fields_idx][0:5] == 'XA:Z:'):
                multi_align_str += line_fields[fields_idx][5:]
                break
        align_list =  multi_align_str[:-1].split(';')
        

    read_seq_len = len(SR_seq)
    for align_str in align_list:
        err_state = False
        fields = align_str.split(',')
        
        ref_seq = LR_seq[fields[0]]
        ref_seq_len = len(ref_seq)
        if (fields[1][0] == '-'):     # Check if seq is reversed complement
            read_seq = SR_seq_rvs_cmplmnt
            pseudo_SR_name = "-" + SR_name
        else:
            read_seq = SR_seq
            pseudo_SR_name = SR_name
        fields[1] = fields[1][1:]
        read_idx = 0
        sub_ref_idx =  1  # 1-offset address
        ref_idx = int(fields[1]) - 1   # convert to 0-offset address
        diff_list = []
        cigar_list = re.split('(M|I|D)', fields[2])
        for idx in range(1, len(cigar_list), 2):
            if (cigar_list[idx - 1].isdigit()):
                if (cigar_list[idx] == 'M'):
                    subseq_len = int(cigar_list[idx - 1])
                    if ((read_idx + subseq_len > read_seq_len) or
                         (ref_idx + subseq_len > ref_seq_len)):
                        err_state = True
                        break
                    read_subseq = numpy.array(list(read_seq[read_idx:(read_idx + subseq_len)]))
                    ref_subseq = numpy.array(list(ref_seq[ref_idx:(ref_idx + subseq_len)]))
                    mut_indices = numpy.where((ref_subseq == read_subseq) == False)[0].tolist()
                    for mut_idx in mut_indices:
                        if (read_subseq[mut_idx] != "N"):
                            diff_list += [str(sub_ref_idx + mut_idx) + ref_subseq[mut_idx] + '>' + read_subseq[mut_idx]]
                    read_idx += subseq_len
                    ref_idx += subseq_len
                    sub_ref_idx += subseq_len
                elif (cigar_list[idx] == 'I'):
                    subseq_len = int(cigar_list[idx - 1])
                    if (read_idx + subseq_len > read_seq_len):
                        err_state = True
                        break
                    insert_str = re.sub(r'N|n', '', read_seq[read_idx:(read_idx + int(cigar_list[idx - 1]))])
                    if (insert_str != ""):
                        diff_list += [str(sub_ref_idx) + '+' + read_seq[read_idx:(read_idx + subseq_len)]]
                    read_idx += subseq_len           
                elif (cigar_list[idx] == 'D'):
                    subseq_len = int(cigar_list[idx - 1])
                    if (ref_idx + subseq_len > ref_seq_len):
                        err_state = True
                        break
                    for del_idx in range(subseq_len):
                        diff_list += [str(sub_ref_idx + del_idx) + '-' + ref_seq[ref_idx + del_idx]]
                    ref_idx += subseq_len
                    sub_ref_idx += subseq_len
            else:
                #print 'Err in cigar: tag : ' + line
                err_state = True
                break
        if ((cigar_list[-1] == '') and (not err_state)):
            if (len(diff_list) == 0):
                diff_list.append("*")
            nav_file.write('\t'.join([pseudo_SR_name, fields[0], fields[1], ' '.join(diff_list),
                                      SR_seq, SR_idx_seq]) + '\n')

sam_file.close()
nav_file.close()
SR_file.close()
SR_idx_file.close()
