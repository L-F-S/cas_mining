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

data_dir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno"

# dataset names
os.chdir("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas")

def find_cas_anno(sample):
    command = "less "+sample+".gff | grep CRISPR >" + data_dir + "/crisprcasanno/"+sample+".crisprcas.gff"
    os.system(command)
    # filter out only proteins
    command2= "less "+ data_dir + "/crisprcasanno/"+sample+".crisprcas.gff | grep protein>" + data_dir + "/justcasanno/"+sample+".cas.gff"
    os.system(command2)
    return

dataset=sys.argv[1]


os.chdir("/scratchCM/tmp_projects/epasolli_darkmatter/allcontigs/"+dataset+"/metabat/genomes_comp50_cont05/prokka/")
print("je suis in "+dataset)
   # tmp=0
for sample in os.listdir():
    if sample.startswith(dataset):
    #        tmp+=1
        os.chdir(sample)
     #       if tmp<4:
      #          print(sample)
        find_cas_anno(sample)
        os.chdir("../")
 
