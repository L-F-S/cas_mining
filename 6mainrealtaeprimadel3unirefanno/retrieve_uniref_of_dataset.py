#!/bin/bash
# 06/06/2019
# Made by L-F-S
# At the University Of Trento, Italy

 # uguale a retrieve_all_annotations.py, ma qui scegli quale
 #dataset passargli (= non cicla su tutti i datasets automaticamente)

 #usage:
     # python retrieve_annotations_of_dataset.py <datasetname>
     # or through sh parall.sh, che lo lancia per tutti

import numpy as np
import pandas as pd
import os
import sys

data_dir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/6mainrealtaeprimadel3unirefanno"

# dataset names
os.chdir("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas")

def find_uniref_cas_anno(binn):
    handle=open(binn)
    from Bio import SeqIO
    for record in SeqIO.parse(handle, "fasta"):
        print("AO")
        if "CRISPR" in record.description:
            line=">"+record.description+"\t"+record.seq+"\t"
        else:
            line=''
        f=open(data_dir+"/"+dataset+"/"+binn+".unirefanno","a")
        print(data_dir+"/"+dataset+"/"+binn+".unirefanno")
        f.write(str(line))
        f.close()
    handle.close()
   # command = "less "+ "binn | grep CRISPR >" + data_dir + "/crisprcasanno/"+binn+".crisprcas.gff"
  #  os.system(command)
    # filter out only proteins
   # command2= "less "+ data_dir + "/crisprcasanno/"+binn+".crisprcas.gff | grep protein>" + data_dir + "/justcasanno/"+binn+".cas.gff"
   # os.system(command2)
    return

dataset=sys.argv[1]

os.chdir("/shares/CIBIO-Storage/CM/scratch/tmp_projects/epasolli_darkmatter/uniref_annotation/"+dataset)
print("je suis in "+dataset)
nbins=0   # tmp=0
for binn in os.listdir():
    if binn.startswith(dataset):
        nbins+=1
        find_uniref_cas_anno(binn) 
print(nbins, "bins evaluated")
#find_uniref_cas_anno("ZeeviD_2015__PNP_DietIntervention_26__bin.10.annotated")
