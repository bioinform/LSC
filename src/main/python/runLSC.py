#!/usr/bin/python

import sys
import os
import datetime
import threading
import string
import commands
from lsccommon import log_command, log_print

################################################################################
def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    path='/'.join(ls[0:-1])+'/'
    if path == "/":
        path = "./"
    return path, filename

def Readcfgfile(cfg_filename):
    results = {}
    cfg = open(cfg_filename,'r')
    for line in cfg:
        line = line.strip()
        if line=='':
            continue
        if not line[0]=='#':
            ls = line.split('=')
            log_print(ls)
            if len(ls)>2:
                log_print('warning: too many = in cfg file')
            results[ls[0].strip()] = ls[1].strip()
    cfg.close()
    return results

################################################################################
if len(sys.argv) >= 2:
    run_pathfilename =  os.path.abspath(os.path.realpath(sys.argv[0]))
    cfg_filename =  sys.argv[1]
else:
    log_print("Correct errors (e.g. homopolymer errors) in long reads, using short read data")
    log_print("usage: python runLSC.py run.cfg")
    log_print("or ./runLSC.py run.cfg")
    sys.exit(1)
################################################################################
version = "1.beta"
python_path = "/usr/bin/python"
mode = 0
Nthread1 = 8
Nthread2 = 8
LR_pathfilename = ''
SR_pathfilename = ''
temp_foldername = 'temp'
output_foldername = 'output'
I_RemoveBothTails = "N"
MinNumberofNonN = "40"
MaxN = "1"
I_nonredundant = "N"
SCD = 20
sort_max_mem = "-1"
clean_up = 0
max_error_rate = "12"
aligner = "bowtie2"
novoalign_options =  "-r All -F FA  -n 300 -o sam" 
bwa_options =  "-n 0.08 -o 10 -e 3 -d 0 -i 0 -M 1 -O 1 -E  1 -N" 
bowtie2_options = "--end-to-end -a -f -L 15 --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.08 --no-unal --omit-sec-seq"
razers3_options = "-i 92 -mr 0 -of sam "

################################################################################
if ((cfg_filename == "-v") or
    (cfg_filename == "--v") or
    (cfg_filename == "-version") or
    (cfg_filename == "--version")):
    log_print("LSC version: " + version)
    exit(0)
    
################################################################################


log_print("=== Welcome to LSC " + version + " ===")

cfg_dt = Readcfgfile(cfg_filename)
for key in cfg_dt:
    if key == "python_path":
        python_path = cfg_dt[key]
    if key == "mode":
        mode = int(cfg_dt[key])
    if key == "Nthread1":
        Nthread1 = int(cfg_dt[key])
    elif key == "Nthread2":
        Nthread2 = int(cfg_dt[key])
    elif key == "LR_pathfilename":
        LR_pathfilename = cfg_dt[key]
    elif key == "LR_filetype":
        LR_filetype = cfg_dt[key]
    elif key == "SR_pathfilename":
        SR_pathfilename = cfg_dt[key]
    elif key == "SR_filetype":
        SR_filetype = cfg_dt[key]
    elif key == "temp_foldername":
        temp_foldername = cfg_dt[key]
    elif key == "output_foldername":
        output_foldername = cfg_dt[key]
    elif key == "I_RemoveBothTails":
        I_RemoveBothTails = cfg_dt[key]
    elif key == "MinNumberofNonN":
        MinNumberofNonN = cfg_dt[key]
    elif key == "MaxN":
        MaxN = cfg_dt[key]
    elif key == "SCD":
        SCD = int(cfg_dt[key])
    elif key == "sort_max_mem":
        sort_max_mem = cfg_dt[key]
    elif key == "I_nonredundant":
        I_nonredundant = cfg_dt[key]
    elif  key ==  "clean_up":
        clean_up =  int(cfg_dt[key])
    elif key == "aligner":
        aligner = cfg_dt[key]
    elif key == "novoalign_options":
        novoalign_options = cfg_dt[key]
    elif key == "bwa_options":
        bwa_options = cfg_dt[key]
    elif key == "bowtie2_options":
        bowtie2_options = cfg_dt[key]
    elif key == "razers3_options":
        razers3_options = cfg_dt[key]
    
################################################################################

if temp_foldername[-1]!='/':
    temp_foldername=temp_foldername+'/'
if output_foldername[-1]!='/':
    output_foldername=output_foldername+'/'

if (not os.path.isdir(output_foldername)):
    log_command('mkdir ' + output_foldername)
if (not os.path.isdir(temp_foldername)):
    if (mode == 2):
        log_print("Error: temp folder does not exist.")
        log_print("Note: You need to run LSC in mode 1 or 0 before running in mode 2.")
        exit(1)
    log_command('mkdir ' + temp_foldername)
if (not os.path.isdir(temp_foldername + "log")):
    log_command('mkdir ' + temp_foldername + "log")


bin_path, run_filename = GetPathAndName(run_pathfilename)
LR_path, LR_filename = GetPathAndName(LR_pathfilename)
SR_path, SR_filename = GetPathAndName(SR_pathfilename)

python_bin_path = python_path + " " + bin_path

t0 = datetime.datetime.now()

################################################################################
if (len(sys.argv) > 2):
    if (sys.argv[2] == "-clean_up"):
        cleanup_cmd = python_bin_path + "clean_up.py " + temp_foldername + " " + str(Nthread1) + " " + str(Nthread2)
        log_command(cleanup_cmd)
    else:
        log_print("")
        log_print("Error: Invalid option " + sys.argv[2])
    exit(0)    
        
################################################################################
# Remove duplicate short reads first
if(SR_filetype == "fa"):
    SR_NL_per_seq = 2
elif(SR_filetype == "fq"):
    SR_NL_per_seq = 4
elif(SR_filetype == "cps"):
    SR_NL_per_seq = -1
else:
    log_print("Err: invalid filetype for short reads")
    exit(1)
    
if ((mode == 0) or 
    (mode == 1)):
    
    if ((I_nonredundant == "N") and
        (SR_filetype != "cps")):
        log_print("=== sort and uniq SR data ===")
        
        fa2seq_cmd = "awk '{if(NR%" + str(SR_NL_per_seq) + "==" + str(2 % SR_NL_per_seq) + ")print $0}' " + SR_pathfilename + " > " + temp_foldername + "SR.seq"
        log_command (fa2seq_cmd)
    
        sort_cmd = "sort -T " + temp_foldername 
        if (sort_max_mem != "-1"):
            sort_cmd += " -S " + str(sort_max_mem)+ " "    
        sort_cmd += " " + temp_foldername + "SR.seq > " + temp_foldername + "SR_sorted.seq"
        log_command(sort_cmd)
    
        uniq_cmd = "uniq -c " + temp_foldername + "SR_sorted.seq > " + temp_foldername + "SR_uniq.seq"
        log_command(uniq_cmd)
    
        uniqseq2fasta_cmd = python_bin_path + "uniqseq2fasta.py " + temp_foldername + "SR_uniq.seq > " + temp_foldername + "SR_uniq.fa"
        log_command(uniqseq2fasta_cmd)
    
        log_print(str(datetime.datetime.now()-t0))
        SR_pathfilename = temp_foldername + "SR_uniq.fa"
        rm_cmd = "rm " + temp_foldername + "SR_uniq.seq " + temp_foldername + "SR_sorted.seq " + temp_foldername + "SR.seq"
        log_command(rm_cmd)

        SR_filetype = "fa"
        
##########################################
# Split the SR FASTQ across CPUs
ext_ls=[]
for i in range(Nthread1):
    ext_ls.append( '.' + string.lowercase[i / 26] + string.lowercase[i % 26] )

if ((mode == 0) or 
    (mode == 1)):
    
        log_print("===split SR:===")    
        if (SR_filetype == "cps"):
            SR_NL = int(commands.getstatusoutput('wc -l ' + SR_pathfilename)[1].split()[0])
            Nsplitline = 1 + (SR_NL / Nthread1)
            if (Nsplitline % 2 == 1):
                Nsplitline +=1
                    
            splitSR_cmd = "split -l " + str(Nsplitline) + " " + SR_pathfilename + " " + temp_foldername + "SR.fa."
            log_command(splitSR_cmd)
            for ext in ext_ls:
                mv_cmd = "mv " + temp_foldername + "SR.fa" + ext + " "  + temp_foldername + "SR.fa" + ext + ".cps"
                log_command(mv_cmd)
            splitSR_cmd = "split -l " + str(Nsplitline/2) + " " + SR_pathfilename.strip()[:-3] + "idx " + temp_foldername + "SR.fa."
            log_command(splitSR_cmd)
            for ext in ext_ls:
                mv_cmd = "mv " + temp_foldername + "SR.fa" + ext + " "  + temp_foldername + "SR.fa" + ext + ".idx"
                log_command(mv_cmd)
                            
            log_print(str(datetime.datetime.now()-t0))
        else:
            SR_NL = int(commands.getstatusoutput('wc -l ' + SR_pathfilename)[1].split()[0])
            Nsplitline = 1 + (SR_NL / Nthread1)
            if ( SR_filetype == "fa"):
                if (Nsplitline % 2 == 1):
                    Nsplitline +=1
            elif ( SR_filetype == "fq"):
                if (Nsplitline % 4 != 0):
                    Nsplitline += (4 - (Nsplitline % 4))
            else:
                log_print("Err: invalid filetype for short reads")
                exit(1)
            
            splitSR_cmd = "split -l " + str(Nsplitline) + " " + SR_pathfilename + " " + temp_foldername + "SR.fa."
            log_command(splitSR_cmd)
                        
            log_print(str(datetime.datetime.now()-t0))

SR_filename = "SR.fa"
if (SR_filetype == "cps"):
    SR_cps_pathfilename = SR_path +  "SR.fa"
else:
    SR_cps_pathfilename = temp_foldername + "SR.fa"

##########################################
# Run HC over the split SR input
if ((mode == 0) or 
    (mode == 1)):

    if (SR_filetype != "cps"):
        log_print("===compress SR.??:===")    
        
        i = 0
        T_compress_SR_ls = []
        for ext in ext_ls:
            compress_SR_cmd = python_bin_path + "compress.py -MinNonN=" + MinNumberofNonN + " -MaxN=" + MaxN + " " + SR_filetype + " " + temp_foldername + SR_filename + ext + " " + temp_foldername + SR_filename + ext + "."
            T_compress_SR_ls.append( threading.Thread(target=log_command, args=(compress_SR_cmd,)) )
            T_compress_SR_ls[i].start()
            i += 1
        for T in T_compress_SR_ls:
            T.join()
        
        log_print(str(datetime.datetime.now()-t0))

        ####################
        # Remove temporary SR split files
        for ext in ext_ls:
            log_command("rm " + temp_foldername + SR_filename + ext)
        ####################

##########################################change output from compress.py and poolchr.py 
# Remove the tails (shorter bits) from the LR, which SHOULD also be the overlapping DNA compliment
if ((mode == 0) or 
    (mode == 1)):

    if I_RemoveBothTails == "Y":   
        log_print("===RemoveBothTails in LR:===")    
        RemoveBothTails_cmd = python_bin_path + "RemoveBothTails.py " + LR_filetype + " " + LR_pathfilename + " " + temp_foldername + "Notwotails_" + LR_filename 
        log_command(RemoveBothTails_cmd)
        LR_filetype = "fa"
        log_print(str(datetime.datetime.now()-t0))

    if I_RemoveBothTails == "Y":
        LR2fa_cmd = python_bin_path + "FASTA2fa.py " + temp_foldername + "Notwotails_" + LR_filename + " " + temp_foldername + "LR.fa"
        deltempLR_cmd = "rm " + temp_foldername + "Notwotails_" + LR_filename  
    else:
        if (LR_filetype == "fa"):
            LR2fa_cmd = python_bin_path + "FASTA2fa.py " + LR_pathfilename + " " + temp_foldername + "LR.fa"
        else:
            LR2fa_cmd = python_bin_path + "FASTQ2fa.py " + LR_pathfilename + " " + temp_foldername + "LR.fa"
            
    log_print(LR2fa_cmd)
    log_command(LR2fa_cmd)
    if I_RemoveBothTails == "Y":
        log_print(deltempLR_cmd)
        log_command(deltempLR_cmd)

LR_filename = "LR.fa"

##########################################
# Compress the long reads
# Build the aligner index and then align short to long 
if ((mode == 0) or 
    (mode == 1)):

    log_print(str(datetime.datetime.now()-t0))   
    
    log_print("===compress LR:===")    
    compress_LR_cmd = python_bin_path + "compress.py -MinNonN=" + MinNumberofNonN + " -MaxN=10000" + " fa " + temp_foldername + LR_filename + " " + temp_foldername + LR_filename +"."
    log_print(compress_LR_cmd)
    log_command(compress_LR_cmd)

    ####################
    LR_NL = int(commands.getstatusoutput('wc -l ' + temp_foldername + "LR.fa.cps")[1].split()[0])
    LR_NR = LR_NL / 2
    ####################
    delLR_cmd = "rm " + temp_foldername + "LR.fa"
    log_print(delLR_cmd)
    log_command(delLR_cmd)
    ####################

    log_print(str(datetime.datetime.now()-t0))

    if (aligner == "bowtie2"):
    
        log_print("===bowtie2 index LR:===")    
        bowtie2_index_cmd = "bowtie2-build -f " + temp_foldername + LR_filename + ".cps " + temp_foldername + LR_filename + ".cps"
        log_command(bowtie2_index_cmd)
        
        log_print(str(datetime.datetime.now()-t0))
        
        
        ##########################################
        log_print("===bowtie2 SR.??.cps:===")    
        
        i=0
        T_bowtie2_ls=[]
        for ext in ext_ls:
            bowtie2_cmd = "bowtie2 " + bowtie2_options + " -x " + temp_foldername + LR_filename + ".cps -U " + temp_foldername + SR_filename + ext + ".cps | samtools view -bS - > " + temp_foldername + SR_filename + ext + ".cps.bam" 
            T_bowtie2_ls.append( threading.Thread(target=log_command, args=(bowtie2_cmd,)) )
            T_bowtie2_ls[i].start()
            i+=1
        for T in T_bowtie2_ls:
            T.join()
    
    elif (aligner == "razers3"):
    
        log_print("===razers3 SR.??.cps:===")    
        
        i=0
        T_razers3_ls=[]
        
        # To be compatible with latest razers3 version 3.1.1: add .fa to *.cps fies.
        razers3_rename_cmd = "mv " + temp_foldername + LR_filename + ".cps " + temp_foldername + LR_filename + ".cps.fa"
        log_command(razers3_rename_cmd)
        for ext in ext_ls:
            # To be compatible with latest razers3 version 3.1.1: add .fa to *.cps fies.
            razers3_rename_cmd = "mv " + temp_foldername + SR_filename + ext + ".cps " + temp_foldername + SR_filename + ext + ".cps.fa"
            log_command(razers3_rename_cmd)
            razers3_cmd = ("razers3 " + razers3_options + " -m " + str(LR_NR) + " -o >(samtools view -bS - >"  + temp_foldername + SR_filename + ext + ".cps.bam) " + 
                           temp_foldername + LR_filename + ".cps.fa " + temp_foldername + SR_filename + ext + ".cps.fa") 
            T_razers3_ls.append( threading.Thread(target=log_command, args=("bash -c \'" + razers3_cmd + "\'",)) )
            T_razers3_ls[i].start()
            i+=1
        for T in T_razers3_ls:
            T.join()
        
        for ext in ext_ls:
            # To be compatible with latest razers3 version 3.1.1: added .fa to *.cps fies.
            razers3_rename_cmd = "mv " + temp_foldername + SR_filename + ext + ".cps.fa " + temp_foldername + SR_filename + ext + ".cps"
            log_command(razers3_rename_cmd)            
        razers3_rename_cmd = "mv " + temp_foldername + LR_filename + ".cps.fa " + temp_foldername + LR_filename + ".cps"
        log_command(razers3_rename_cmd)
        
    elif (aligner == "bwa"):
    
        log_print("===bwa index LR:===")    
        bwa_index_cmd = "bwa index " + temp_foldername + LR_filename + ".cps"
        log_command(bwa_index_cmd)
        
        log_print(str(datetime.datetime.now()-t0))
        
        
        ##########################################
        log_print("===bwa aln SR.??.cps:===")    
        
        i=0
        T_bwa_ls=[]
        for ext in ext_ls:
            bwa_cmd = "bwa aln " + bwa_options + " " + temp_foldername + LR_filename + ".cps " + temp_foldername + SR_filename + ext + ".cps > " + temp_foldername + SR_filename + ext + ".cps.sai" 
            T_bwa_ls.append( threading.Thread(target=log_command, args=(bwa_cmd,)) )
            T_bwa_ls[i].start()
            i+=1
        for T in T_bwa_ls:
            T.join()
            
        ##########################################
        log_print("===bwa samse SR.??.cps.bam :===")    
        
        i=0
        T_bwa_ls=[]
        for ext in ext_ls:
            bwa_cmd = "bwa samse -n " + str(LR_NR) + " " + temp_foldername + LR_filename + ".cps " + temp_foldername + SR_filename + ext + ".cps.sai " + temp_foldername + SR_filename + ext + ".cps | samtools view -bS - > " + temp_foldername + SR_filename + ext + ".cps.bam" 
            T_bwa_ls.append( threading.Thread(target=log_command, args=(bwa_cmd,)) )
            T_bwa_ls[i].start()
            i+=1
        for T in T_bwa_ls:
            T.join()
            
        ####################
        for ext in ext_ls:
            log_command("rm " + temp_foldername + SR_filename + ext + ".cps.sai")
        ####################
    
    else:
    
        log_print("===novoindex LR:===")    
        novoindex_cmd = "novoindex " + temp_foldername + LR_filename + ".cps.nix " + temp_foldername + LR_filename + ".cps"
        log_command(novoindex_cmd)
        
        log_print(str(datetime.datetime.now()-t0))
        
        
        ##########################################
        log_print("===novoalign SR.??.cps:===")    
        
        i=0
        T_novoalign_ls=[]
        for ext in ext_ls:
            novoalign_cmd = "novoalign " + novoalign_options + " -r Ex " + str(LR_NR) + " -d " + temp_foldername + LR_filename + ".cps.nix -f " + temp_foldername + SR_filename + ext + ".cps > " + temp_foldername + SR_filename + ext + ".cps.sam" 
            T_novoalign_ls.append( threading.Thread(target=log_command, args=(novoalign_cmd,)) )
            T_novoalign_ls[i].start()
            i+=1
        for T in T_novoalign_ls:
            T.join()
    
            
    log_print(str(datetime.datetime.now()-t0))

    ##########################################
    # Convert the SAM file to a NAV file
    log_print("===samParser SR.??.cps.nav:===")
    
    i=0
    T_samParser_ls=[]
    for ext in ext_ls:
        samParser_cmd = (python_bin_path + "samParser.py " + temp_foldername + LR_filename + ".cps " + temp_foldername + SR_filename + ext +
                         " <(samtools view -h " + temp_foldername + SR_filename + ext + ".cps.bam) " + max_error_rate + " ") 
        if (aligner == "bwa"):
            samParser_cmd += " F "   # Setting one_line_per_alignment parameter 
        else:
            samParser_cmd += " T " 
        samParser_cmd += " 2> " + temp_foldername + SR_filename + ext + ".cps.samParser.log"
        samParser_cmd += " | gzip -1 > " + temp_foldername + SR_filename + ext + ".cps.nav.gz"
    
        T_samParser_ls.append( threading.Thread(target=log_command, args=("bash -c \'" + samParser_cmd + "\'",)) )
        T_samParser_ls[i].start()
        i+=1
    i = 0
    for T in T_samParser_ls:
        T.join()

            
    log_print(str(datetime.datetime.now()-t0))
    
##########################################
# Build complete CPS and IDX files
    if (SR_filetype != "cps"): 
        if (clean_up < 2):
            log_print("===cat SR.??.cps:===")    
            for ext in ext_ls:
                log_command( "cat " + temp_foldername + SR_filename + ext + ".cps >> " + temp_foldername + SR_filename + ".cps", printcommand=False )
                log_command( "rm " + temp_foldername + SR_filename + ext + ".cps ", printcommand=False )
            log_print(str(datetime.datetime.now()-t0))
            ####################
            log_print("===cat SR.??.idx:===")    
            for ext in ext_ls:
                log_command( "cat " + temp_foldername + SR_filename + ext + ".idx >> " + temp_foldername + SR_filename + ".idx", printcommand=False )
                log_command( "rm " + temp_foldername + SR_filename + ext + ".idx ", printcommand=False )
            log_print(str(datetime.datetime.now()-t0))
        elif (clean_up == 2):
            for ext in ext_ls:
                log_command("rm " + temp_foldername + SR_filename + ext + ".idx")
                log_command("rm " + temp_foldername + SR_filename + ext + ".cps")
    else:
        for ext in ext_ls:
            log_command("rm " + temp_foldername + SR_filename + ext + ".cps")
            log_command("rm " + temp_foldername + SR_filename + ext + ".idx")
        ####################

    # Build complete SAM and NAV files
    ####################
    if (clean_up < 2):
        log_print("===cat SR.??.cps.bam :===")    
        log_command("touch " + temp_foldername + SR_filename + ".cps.bam")
        for ext in ext_ls:
            log_command("cat " + temp_foldername + SR_filename + ext + ".cps.bam" + " >> " + temp_foldername + SR_filename + ".cps.bam", printcommand=False )
            log_command("rm " + temp_foldername + SR_filename + ext + ".cps.bam")
        log_print(str(datetime.datetime.now()-t0))
    elif (clean_up == 2):
        # Remove bam, alignment summary files per thread
        for ext in ext_ls:
            log_command("rm " + temp_foldername + SR_filename + ext + ".cps.bam")
    
    ####################
    # For uncompressed version
    #log_print("===cat SR.??.cps.nav :===")
    #log_command("touch " + temp_foldername + SR_filename + ".cps.nav")
    #for ext in ext_ls:
        #log_command("cat " + temp_foldername + SR_filename + ext + ".cps.nav" + " >> " + temp_foldername + SR_filename + ".cps.nav", printcommand=False )
        #log_command("rm " + temp_foldername + SR_filename + ext + ".cps.nav")
    #log_print(str(datetime.datetime.now()-t0))
    # For uncompressed version
    log_print("===cat SR.??.cps.nav :===")
    nav_filename_ls = []
    for ext in ext_ls:
        nav_filename_ls.append(temp_foldername + SR_filename + ext + ".cps.nav.gz")
    log_command("cat " + " ".join(nav_filename_ls) + " > " + temp_foldername + SR_filename + ".cps.nav.gz")
    log_command("rm " + " ".join(nav_filename_ls))
    log_print(str(datetime.datetime.now()-t0))
    ####################
    for ext in ext_ls:
        log_command("mv " + temp_foldername + "SR.fa" + ext + ".cps.samParser.log " + temp_foldername + "log")
    ####################

####################

# Return after alignment in case of mode 1
if (mode == 1):
    exit(0)


##########################################
log_print("===genLR_SRmapping SR.??.cps.nav:===")    
    
genLR_SRmapping_cmd = python_bin_path + "genLR_SRmapping.py "  + temp_foldername + " " +  temp_foldername + SR_filename + ".cps.nav.gz " + temp_foldername + LR_filename
genLR_SRmapping_cmd += " " + str(SCD) + " " + str(Nthread2) + " " + str(sort_max_mem) 
log_command(genLR_SRmapping_cmd)
if (clean_up == 2):
    log_command("rm " + temp_foldername + SR_filename + ".cps.nav.gz")
log_print(str(datetime.datetime.now()-t0))

##########################################
log_print("===split LR_SR.map:===")    

LR_SR_map_NR = int(commands.getstatusoutput('wc -l ' + temp_foldername +"LR_SR.map")[1].split()[0])

if (LR_SR_map_NR == 0):
    log_print("Error: No short reads was aligned to long read. LSC could not correct any long read sequence.")
    exit(1)
    
Nsplitline = 1 + (LR_SR_map_NR/Nthread2)

Nthread2_temp = int(LR_SR_map_NR)/int(Nsplitline)
if ((LR_SR_map_NR % Nsplitline) != 0):
    Nthread2_temp += 1
if (Nthread2_temp < Nthread2):
    Nthread2 = Nthread2_temp

ext2_ls=[]
for i in range(Nthread2):
    ext2_ls.append( '.' + string.lowercase[i / 26] + string.lowercase[i % 26] )
    
splitLR_SR_map_cmd = "split -l " + str(Nsplitline) + " " + temp_foldername + "LR_SR.map" + ' ' + temp_foldername + "LR_SR.map" +"."
log_command(splitLR_SR_map_cmd)
if (clean_up == 2):
    log_command("rm " + temp_foldername + "LR_SR.map ")

log_print(str(datetime.datetime.now()-t0))

##########################################
log_print("===correct.py LR_SR.map.??_tmp :===")    

log_print(str(datetime.datetime.now()-t0))

i=0
T_correct_for_piece_ls=[]
for ext in ext2_ls:
    correct_for_piece_cmd = python_bin_path + "correct_nonredundant.py " + temp_foldername + "LR_SR.map" + ext  + " " + temp_foldername + 'LR.fa.readname  > ' + temp_foldername + "LR_SR.map_emtry_ls" + ext
    T_correct_for_piece_ls.append( threading.Thread(target=log_command, args=(correct_for_piece_cmd,)) )
    T_correct_for_piece_ls[i].start()
    i+=1
for T in T_correct_for_piece_ls:
    T.join()

log_print(str(datetime.datetime.now()-t0))

####################
if (clean_up >= 1):
    for ext in ext2_ls:
        delLR_SR_map_aa_tmp_cmd = "rm " + temp_foldername + "LR_SR.map" + ext 
        log_command(delLR_SR_map_aa_tmp_cmd)
####################

##########################################

log_print("===cat full_LR_SR.map.fa :===")    

temp_filename_ls = []
for ext in ext2_ls:
    temp_filename_ls.append( temp_foldername + "full_LR_SR.map" + ext )
log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldername + "full_LR.fa" )

log_print("===cat corrected_LR_SR.map.fa :===")    

temp_filename_ls = []
for ext in ext2_ls:
    temp_filename_ls.append( temp_foldername + "corrected_LR_SR.map" + ext )
log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldername + "corrected_LR.fa" )

log_print("===cat corrected_LR_SR.map.fq :===")    

temp_filename_ls = []
for ext in ext2_ls:
    temp_filename_ls.append( temp_foldername + "corrected_LR_SR.map" + ext + '.fq' )
log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldername + "corrected_LR.fq" )

log_print("===cat uncorrected_LR_SR.map.fa :===")    

temp_filename_ls = []
for ext in ext2_ls:
    temp_filename_ls.append( temp_foldername + "uncorrected_LR_SR.map" + ext  )
log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldername + "uncorrected_LR.fa" )

####################
if (clean_up >= 1):
    for ext in ext2_ls:
        del_LR_SR_coverage_aa_cmd = "rm " + temp_foldername + "corrected_LR_SR.map" + ext + ".fq"
        log_command(del_LR_SR_coverage_aa_cmd)
        delfull_LR_SR_map_aa_fa_cmd = "rm " + temp_foldername + "full_LR_SR.map" + ext 
        log_command(delfull_LR_SR_map_aa_fa_cmd)
        delcorr_LR_SR_map_aa_fa_cmd = "rm " + temp_foldername + "corrected_LR_SR.map" + ext
        log_command(delcorr_LR_SR_map_aa_fa_cmd)
        deluncorr_LR_SR_map_aa_fa_cmd = "rm " + temp_foldername + "uncorrected_LR_SR.map" + ext
        log_command(deluncorr_LR_SR_map_aa_fa_cmd)
####################

####################
for ext in ext2_ls:
    log_command("mv " + temp_foldername + "LR_SR.map_emtry_ls" + ext + " " + temp_foldername + "log")
####################


