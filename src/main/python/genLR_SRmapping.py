#!/usr/bin/python

import sys
import string
import commands
import threading
import random
import gzip
from lsccommon import log_command, log_print


################################################################################
# Debug flags
printSCD = False

################################################################################

if len(sys.argv) >= 2:
    temp_foldername = sys.argv[1]
    nav_filename = sys.argv[2]
    LR_filename = sys.argv[3]
    SR_cvrg_threshold = int(sys.argv[4])
    Nthread = int(sys.argv[5])
    sort_max_mem = sys.argv[6]
    
    
else:
    log_print("usage: python convertNAV.py temp_foldername LR_filename nav_filename Nthread sort_max_mem")
    log_print("or ./convertNAV.py temp_foldername LR_filename nav_filename Nthread sort_max_mem")
    sys.exit(1)

################################################################################
# Splitting the nav file
sort_cmd = "bash -c 'sort -T " + temp_foldername
if (sort_max_mem != "-1"):
    sort_cmd += " -S " + str(sort_max_mem) + " "
sort_cmd += " -nk 2 --compress-program=bzip2 --parallel=" + str(Nthread) + " <(zcat "  + nav_filename + " ) | gzip -1 > " + nav_filename + ".sort.gz'"

log_command(sort_cmd)

log_print("Done with sorting")

################################################################################
LR_cps_file = open(LR_filename + '.cps','r')
LR_cps_dict={}
for line in LR_cps_file:
    if line[0]=='>':
        readname = line[1:-1]
    else:
        LR_cps_dict[readname] = line.strip()
LR_cps_file.close()

LR_idx_file = open(LR_filename +'.idx','r')
LR_idx_dict={}
for line in LR_idx_file:
    fields=line.strip().split('\t')
    LR_idx_dict[fields[0]] = '\t'.join(fields[1:])
LR_idx_file.close()

LR_name_file = open(LR_filename +'.readname','r')
LR_name_dict={}
for line in LR_name_file:
    fields=line.strip().split('\t')
    LR_name_dict[fields[0]] = fields[1]
LR_name_file.close()

nav_file= gzip.open(nav_filename + ".sort.gz" ,'r')
nav_cvrg_file= gzip.open(nav_filename + ".sort.gz" ,'r')   # This used to compute coverage
LR_SR_mapping_file = open(temp_foldername + "LR_SR.map",'w')
if (printSCD):
    LR_SR_coverage_file = open(temp_foldername + "LR_SR.scd",'w')
    LR_uSR_coverage_file = open(temp_foldername + "LR_SR.uscd",'w')
    LR_SR_coverage_selected_file = open(temp_foldername + "LR_SR.scd.selected",'w')
    LR_uSR_coverage_selected_file = open(temp_foldername + "LR_SR.uscd.selected",'w')

def write_2_LR_SR_map_file(nav_file, num_lines, LR_name,
                           LR_coverage_list, LR_uniq_coverage_list):
    
    # Store LR-SR SCD
    if (printSCD):
        LR_SR_coverage_file.write(">" + LR_name_dict[LR_name] + "\n")
        LR_coverage_str_list = [str(i) for i in LR_coverage_list]
        LR_SR_coverage_file.write(",".join(LR_coverage_str_list) + "\n")
        
        LR_uSR_coverage_file.write(">" + LR_name_dict[LR_name] + "\n")
        LR_coverage_str_list = [str(i) for i in LR_uniq_coverage_list]
        LR_uSR_coverage_file.write(",".join(LR_coverage_str_list) + "\n")
        
        LR_coverage_str_list_temp = [0] * len(LR_coverage_str_list)
        LR_uniq_coverage_str_list_temp = [0] * len(LR_coverage_str_list)
    
    LR_SR_list = []
    
    for line_number in range(num_lines):
        line = nav_file.readline()
        if line[0]=='#':
            continue
        line_list=line.strip().split('\t')
        
        if (LR_name != line_list[1] ):
            print "Error: Unexpected LR name: " + line_list[1]
            exit(1)
        SR_rpt = int(line_list[0].split('_')[1])
        pos = int(line_list[2])
        SR_len = len(line_list[3])
        region_cvrg = min(LR_coverage_list[pos:(pos+SR_len)])
        if (SR_cvrg_threshold < 0):    # Disable the feature for negative threshold
            SR_prob = 1.1
        else:
            SR_prob = 1./region_cvrg * SR_cvrg_threshold   # Note: SR_prob could be greater than 1
        
        use_SR = False
        for SR_rpt_idx in range(SR_rpt):
            if (random.random() <= SR_prob):
                use_SR = True
                break
        if (use_SR):
            if (printSCD):
                LR_coverage_str_list_temp[pos:(pos+SR_len)] = [i + SR_rpt for i in  LR_coverage_str_list_temp[pos:(pos+SR_len)]]
                LR_uniq_coverage_str_list_temp[pos:(pos+SR_len)] = [i + 1 for i in  LR_uniq_coverage_str_list_temp[pos:(pos+SR_len)]]
            if (line_list[3] == "*"):
                line_list[3] = ""
            line_list_temp = [line_list[0],line_list[2]] + line_list[3:]
            LR_SR_list.append(line_list_temp)
    if (printSCD):
        LR_SR_coverage_selected_file.write(">" + LR_name_dict[LR_name] + "\n")
        LR_coverage_str_list = [str(i) for i in LR_coverage_str_list_temp]
        LR_SR_coverage_selected_file.write(",".join(LR_coverage_str_list) + "\n")
        
        LR_uSR_coverage_selected_file.write(">" + LR_name_dict[LR_name] + "\n")
        LR_coverage_str_list = [str(i) for i in LR_uniq_coverage_str_list_temp]
        LR_uSR_coverage_selected_file.write(",".join(LR_coverage_str_list) + "\n")

    temp_SR_ls = []
    ls_SR_seq = []
    ls_SR_idx_seq = []
    for SR in LR_SR_list:
        temp_SR_ls.append(SR[0]+','+ SR[1]+','+SR[2])
        ls_SR_seq.append(SR[3])
        if (len(SR) > 4):
            ls_SR_idx_seq.append('\t'.join([SR[4], SR[5]]))
        else:
            ls_SR_idx_seq.append('\t'.join(["", ""]))   # no compression point

    input_ls= [LR_cps_dict[prev_LR_name], LR_idx_dict[prev_LR_name], ';'.join(temp_SR_ls),prev_LR_name,'kinfai'.join(ls_SR_seq), 'kinfai'.join(ls_SR_idx_seq)]
    LR_SR_mapping_file.write('yue'.join(input_ls)+'\n')


line_list = nav_cvrg_file.readline().strip().split('\t')
LR_name = line_list[1]
LR_coverage_list  = [0] * len(LR_cps_dict[LR_name])
LR_uniq_coverage_list = [0] * len(LR_cps_dict[LR_name])
SR_rpt = int(line_list[0].split('_')[1])
pos = int(line_list[2])
SR_len = len(line_list[3])
LR_coverage_list[pos:(pos+SR_len)] = [(i + SR_rpt) for i in LR_coverage_list[pos:(pos+SR_len)]]
LR_uniq_coverage_list[pos:(pos+SR_len)] = [(i + 1) for i in LR_uniq_coverage_list[pos:(pos+SR_len)]]
prev_LR_name = LR_name
num_lines = 0   # pre-incremented

for line in nav_cvrg_file:
    num_lines += 1
    if line[0]!='#':
         
        line_list=line.strip().split('\t')

        LR_name = line_list[1]
        if (prev_LR_name != LR_name):
            write_2_LR_SR_map_file(nav_file, num_lines, prev_LR_name, LR_coverage_list, LR_uniq_coverage_list)
            
            LR_coverage_list  = [0] * len(LR_cps_dict[LR_name])
            LR_uniq_coverage_list  = [0] * len(LR_cps_dict[LR_name])
            num_lines = 0
            prev_LR_name = LR_name

        SR_rpt = int(line_list[0].split('_')[1])
        pos = int(line_list[2])
        SR_len = len(line_list[3])
        LR_coverage_list[pos:(pos+SR_len)] = [(i + SR_rpt) for i in LR_coverage_list[pos:(pos+SR_len)]]
        LR_uniq_coverage_list[pos:(pos+SR_len)] = [(i + 1) for i in LR_uniq_coverage_list[pos:(pos+SR_len)]]

num_lines += 1
write_2_LR_SR_map_file(nav_file, num_lines, prev_LR_name, LR_coverage_list, LR_uniq_coverage_list)
    
nav_cvrg_file.close()
nav_file.close()
LR_SR_mapping_file.close()
if (printSCD):
    LR_SR_coverage_file.close()
    LR_uSR_coverage_file.close()
    LR_SR_coverage_selected_file.close()
    LR_uSR_coverage_selected_file.close()
    
delSRnavsort_cmd = "rm " + nav_filename + ".sort.gz"
log_command(delSRnavsort_cmd)

log_print("Done with generating LR_SR.map file")

