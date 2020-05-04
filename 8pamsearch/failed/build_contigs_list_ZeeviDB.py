#!/bin/bash
# Sun Jul 28 10:40:55 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy
# builds a comma separated file listing all contigs, and, for every contig, the Genome Name, Study, Sample Name, SGB ID 
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, '/home/lorenzo.signorini/utils/')
import filename_discrepancies


##################### MAIN #################################################
# sbagliato: piu veloce se lo faccio al contrario., bin by bin.
dataset="ZeeviD_2015_B" 
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" 
originaldir="/scratchCM/tmp_projects/epasolli_darkmatter/allcontigs/"+dataset+"/metabat/genomes_comp50_cont05/prokka"
S3_alternative_dataset_name="ZeeviD_2015"

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
       gff_file=open(binn+".gff","r")
       for line in gff_file.readlines():
           if line.startswith("##sequence"):
               contigname=line.strip("\n").split(" ")[1]
               binned_contigs_summary_for_dataset.loc[index]=[contigname,genomewithoutmegahit,S3_alternative_dataset_name,justsample, int(SGB_ID)]
               index+=1

               # keep track of binned contigs for a later check of unbinned conts
               binned_contigs_of_sample[justsample][0].append(contigname)
               if not line.startswith("#"):
                   break
       gff_file.close()
    
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
    f=open(allcontsdir+"/"+samplename_originalfile+"/contigs_filtered.fasta","r")
    for line in f.readlines():
        if line.startswith(">"):
            list_of_all_contigs_for_that_sample.append(line.lstrip(">").rstrip("\n"))

    f.close()
    unbinned_contigs=set(list_of_all_contigs_for_that_sample)-set(binned_contigs_of_sample[justsample][0])
    for contig in unbinned_contigs:
        unbinned_contigs_summary_for_dataset.loc[uind]=[contig,S3_alternative_dataset_name+"__"+justsample+"__unbinned",S3_alternative_dataset_name,justsample,0]
        uind+=1

# 3.merge datasets and  write to  file
contigs_summary_for_dataset=binned_contigs_summary_for_dataset.append(unbinned_contigs_summary_for_dataset, ignore_index=True)
print("writing to file")
contigs_summary_for_dataset.to_csv(datadir+"/4exploratoryanalyses/pamsearch/contigs_summary_"+dataset+".csv")
