#!/usr/bin/python

import sys
import os
import commands
import string

################################################################################
def log_print(print_str):
    os.system("echo " + str(print_str))

def rm_command(filename):
    if (os.path.isfile(filename)):
        os.system("rm " + filename)
        
################################################################################
if len(sys.argv) >= 4:
    temp_foldername =  sys.argv[1]
    Nthread1 =  int(sys.argv[2])
    Nthread2 =  int(sys.argv[3])
else:
    log_print("usage: python clean_up.py temp_foldername Nthread1 Nthread2 ")
    log_print("or ./clean_up.py  temp_foldername Nthread1 Nthread2")
    sys.exit(1)

ext_ls=[]
for i in range(Nthread1):
    ext_ls.append( '.' + string.lowercase[i / 26] + string.lowercase[i % 26] )
    
SR_filename = "SR.fa"
for ext in ext_ls:
    rm_command(temp_foldername + SR_filename + ext + ".cps.nav")
    rm_command(temp_foldername + SR_filename + ext + ".cps.sam")
    rm_command(temp_foldername + SR_filename + ext + ".idx")
    rm_command(temp_foldername + SR_filename + ext + ".cps")

ext2_ls=[]
for i in range(Nthread2):
    ext2_ls.append( '.' + string.lowercase[i / 26] + string.lowercase[i % 26] )


for ext in ext2_ls:
    rm_command(temp_foldername + "LR_SR.map" + ext)

    rm_command(temp_foldername + "corrected_LR_SR.map" + ext + ".fq")
    rm_command(temp_foldername + "full_LR_SR.map" + ext )
    rm_command(temp_foldername + "corrected_LR_SR.map" + ext )
    rm_command(temp_foldername + "uncorrected_LR_SR.map" + ext)
####################
