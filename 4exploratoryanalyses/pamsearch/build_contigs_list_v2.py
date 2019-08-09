#!/bin/bash
# Sun Jul 28 10:40:55 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy
# builds a comma separated file listing all contigs, and, for every contig, the Genome Name, Study, Sample Name, SGB ID 
# added subprocess module to handle fileopening with grep

import os
import sys
import numpy as np
import pandas as pd
import subprocess

sys.path.insert(0, '/home/lorenzo.signorini/utils/')
import filename_discrepancies


##################### MAIN #################################################
dataset=sys.argv[1] 
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" 
originaldir="/scratchCM/tmp_projects/epasolli_darkmatter/allcontigs/"+dataset+"/metabat/genomes_comp50_cont05/prokka"
S3_alternative_dataset_name=filename_discrepancies.s3(dataset)

print("processing dataset "+ dataset)
print("S3 datasetname = ", S3_alternative_dataset_name)
s3=pd.read_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/S3Segata.csv",index_col=0 )

allcontsdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/epasolli_metagenomicassembly/"+"/"+dataset

# 1. retrieve binned contigs
print("+++++++++++++++retrieving binned contigs")
binned_contigs_summary_for_dataset=pd.DataFrame(columns=["Contig Name", "Genome Name","Study","Sample Name","SGB ID"])
index=0
binned_contigs_of_sample={}
os.chdir(originaldir)
for binn in os.listdir():
    print("binn ",binn)
    if not binn.startswith("summary"):
        os.chdir(originaldir+"/"+binn)
        
        justsample= binn.split("__")[1].rstrip("_megahit")
        genomewithoutmegahit=S3_alternative_dataset_name+"__"+justsample+"__"+binn.split("__")[2]
        SGB_ID = str(s3[s3["Genome Name"]==genomewithoutmegahit]["SGB ID"]).split()[1]
 
        if not justsample in binned_contigs_of_sample.keys():
            binned_contigs_of_sample[justsample]=[[],0]
        if "megahit" in binn:
            binned_contigs_of_sample[justsample][1]=True
        else: 
            binned_contigs_of_sample[justsample][1]=False
        
        gff_file=binn+".gff"
        gff_cmd="grep '##sequence' "+gff_file 
        contigs_unpolished=subprocess.check_output(gff_cmd, shell=True, universal_newlines=True).strip("\n")
        for unpolished_contig in contigs_unpolished.split("\n"):
            contigname=unpolished_contig.lstrip("##sequence-region ")
            binned_contigs_summary_for_dataset.loc[index]=[contigname,genomewithoutmegahit,S3_alternative_dataset_name,justsample, int(SGB_ID)]
            index+=1
 
            # keep track of binned contigs for a later check of unbinned conts
            binned_contigs_of_sample[justsample][0].append(contigname)

    
# 2. retrieve unbinned contigs
unbinned_contigs_summary_for_dataset=pd.DataFrame(columns=["Contig Name", "Genome Name","Study","Sample Name","SGB ID"])
uind=0
print("+++++++++++++++++++++++++++++++retrieving unbinned contigs")
for justsample in binned_contigs_of_sample.keys():
    print("sample", justsample)
    if binned_contigs_of_sample[justsample][1]:
        samplename_originalfile=justsample+"_megahit"
    else:
        samplename_originalfile=justsample
    list_of_all_contigs_for_that_sample=[]
    filename=allcontsdir+"/"+samplename_originalfile+"/contigs_filtered.fasta"
    cmd="grep '^>' "+filename
    contigs_unpolished=subprocess.check_output(cmd, shell=True, universal_newlines=True).strip("\n")
    for unpolished_contig in contigs_unpolished.split("\n"):
        contigname=unpolished_contig.lstrip(">")
        list_of_all_contigs_for_that_sample.append(contigname)
    unbinned_contigs=set(list_of_all_contigs_for_that_sample)-set(binned_contigs_of_sample[justsample][0])
    for contig in unbinned_contigs:
        unbinned_contigs_summary_for_dataset.loc[uind]=[contig,S3_alternative_dataset_name+"__"+justsample+"__unbinned",S3_alternative_dataset_name,justsample,0]
        uind+=1

# 3.merge datasets and  write to  file
contigs_summary_for_dataset=binned_contigs_summary_for_dataset.append(unbinned_contigs_summary_for_dataset, ignore_index=True)
print("writing to file")
contigs_summary_for_dataset.to_csv(datadir+"/4exploratoryanalyses/pamsearch/contigs_summary_"+dataset+".csv")
