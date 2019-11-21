# Thu Nov 21 14:10:26 CET 2019
# Made by L-F-S
# At the University Of Trento, Italy

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies


dataset="SmitsSA_2017" #sys.argv[1] #TODO switch test dataset
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" #TODO add nstepdir
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"


def unknownSGB(data,UI):
    n=0
    while not (UI.lower()=="u" or UI.lower()=="k" or UI.lower()=="b"):
        UI=str(input("Invalid option. Please choose one between [u,k,b] for unknown, known, both"))
        n+=1
        if n>9:
            r=np.random.randint(1,10)
            if r<5:
                print("eddaje n po'. u k b EBBASTA. Va bene pure maiuscolo non è difficile su")
                n=0
    if UI=="u":
        return data[data.uSGB=="Yes"]
    elif UI=="k":
        return data[data.uSGB=="No"]
    else:
        return data

