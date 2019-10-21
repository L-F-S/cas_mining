#!/bin/bash
# Fri Aug 2 09:59:19 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy
# Usage:

# python 1build_known_effector_cas_table.py <datasetname> <effector_cas_name>
# <effector_cas_name> can be Cas9 or Cpf1 (as of 17/09/2019)

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies


wdataset=sys.argv[1] #TODO switch test dataset  use folderloop, NOT folderloops  use folderloop, NOT folderloops33
s3dataset=filename_discrepancies.s3(wdataset)
tabledir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/3tabellazza"
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" #TODO add nstepdir
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"

DF=pd.read_csv(tabledir+"/"+ wdataset +"/crisprcas_hits_table_"+wdataset+".csv", index_col=0)
#set +"/crisprcasExtractct sequence information from GENOMES WITH <cas9> annotation.
feature=sys.argv[2]  #TOD  #TODOO
#SGBs=pd.read_excel(datadir+"S4Segata.xlsx")
# OPTIONAL TODO : subsample for a single SGB
#chosen_SGB=15286  # Cibiobacter SGB
#single_SGB=tabellazza[tabellazza["SGB ID"]==chosen_SGB]
#subset_of_genomes=single_SGB 

subset_of_genomes=DF

genomes_with_that_feature=subset_of_genomes[subset_of_genomes.prokka_cas.str.contains(feature)==True]
print("Found ", genomes_with_that_feature.shape[0], " genomes with that feature.")
genomes_with_that_feature.to_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/"+s3dataset+"/"+s3dataset+"_known_"+feature+"_genomes.csv")
##############################

