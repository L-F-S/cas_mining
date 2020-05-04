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
dataset=sys.argv[1] 
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" 
originaldir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/epasolli_darkmatter/allcontigs/"+dataset+"/metabat/genomes_comp50_cont05/prokka"
S3_alternative_dataset_name=filename_discrepancies.s3(dataset)

print("processing dataset "+ dataset)
print("S3 datasetname = ", S3_alternative_dataset_name)
s3=pd.read_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/S3Segata_split/S3_subs_"+S3_alternative_dataset_name+".csv",index_col=0 )

allcontsdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/epasolli_metagenomicassembly/"+dataset

# 1. retrieve binned contigs
print("+++++++++++++++retrieving binned contigs")
binned_contigs_summary_for_dataset=pd.DataFrame(columns=["Contig Name", "Genome Name","Study","Sample Name","SGB ID"],index=range(s3.shape[0]))
binned_contigs_of_sample={}
os.chdir(originaldir)
index=0
for binn in os.listdir():
   if not binn.startswith("summary"):
       print("binn ",binn)
#       os.chdir(originaldir+"/"+binn)
       
       justsample= binn.split("__")[1].rstrip("_megahit")
       genomewithoutmegahit=S3_alternative_dataset_name+"__"+justsample+"__"+binn.split("__")[2]
       SGB_ID = str(s3[s3["Genome Name"]==genomewithoutmegahit]["SGB ID"]).split()[1]

       if not justsample in binned_contigs_of_sample.keys():
           binned_contigs_of_sample[justsample]=[[],0]
       if "megahit" in binn:
           binned_contigs_of_sample[justsample][1]=True
       else: 
           binned_contigs_of_sample[justsample][1]=False
       gff_file=open(originaldir+"/"+binn+"/"+binn+".gff","r")
       for line in gff_file.readlines():
           if line.startswith("##sequence"):
               contigname=line.strip("\n").split(" ")[1]
               binned_contigs_summary_for_dataset.set_value(index,"Contig Name",contigname)
               binned_contigs_summary_for_dataset.set_value(index,"Genome Name",genomewithoutmegahit)
               binned_contigs_summary_for_dataset.set_value(index,"Study",S3_alternative_dataset_name)
               binned_contigs_summary_for_dataset.set_value(index,"Sample Name",justsample)
               binned_contigs_summary_for_dataset.set_value(index,"SGB ID",SGB_ID)
               index+=1
               # keep track of binned contigs for a later check of unbinned conts
               binned_contigs_of_sample[justsample][0].append(contigname)
           if not line.startswith("#"):
               break
       gff_file.close()
    
# 2. retrieve unbinned contigs
print("+++++++++++++++++++++++++++++++retrieving unbinned contigs")
firstsample=True
for justsample in binned_contigs_of_sample.keys():
    print("sample", justsample)
    if binned_contigs_of_sample[justsample][1]:
        samplename_originalfile=justsample+"_megahit"
    else:
        samplename_originalfile=justsample
    list_of_all_contigs_for_that_sample=[]
    print("gathering all contigs")
    f=open(allcontsdir+"/"+samplename_originalfile+"/contigs_filtered.fasta","r")
    for line in f.readlines():
        if line.startswith(">"):
            list_of_all_contigs_for_that_sample.append(line.lstrip(">").rstrip("\n"))

    f.close()
    print("subtracting binned contigs")
    unbinned_contigs=set(list_of_all_contigs_for_that_sample)-set(binned_contigs_of_sample[justsample][0])
    unbinned_contigs_summary_for_dataset_for_sample=pd.DataFrame(columns=["Contig Name", "Genome Name","Study","Sample Name","SGB ID"] index=range(len(unbinned_contigs)))
    print("creating dataset of unbinned contigs")
    for uind,  contig in enumerate(unbinned_contigs):
        unbinned_contigs_summary_for_dataset_for_sample[uind]=[contig,S3_alternative_dataset_name+"__"+justsample+"__unbinned",S3_alternative_dataset_name,justsample,0]

    # append dataset from unbinned samples
    print("appending data")
    if firstsample
        unbinned_contigs_summary_for_dataset=unbinned_contigs_summary_for_dataset_for_sample
        firstsample=False
    else
        unbinned_contigs_summary_for_dataset.append(unbinned_contigs_summary_for_dataset_for_sample , ignore_index=True)

# 3.merge datasets and  write to  file
print("appending binned + unbinned datasets")
contigs_summary_for_dataset=binned_contigs_summary_for_dataset.append(unbinned_contigs_summary_for_dataset, ignore_index=True)
print("writing to file")
contigs_summary_for_dataset.to_csv(datadir+"/4exploratoryanalyses/pamsearch/contigs_summary_"+dataset+".csv")
