# Thu Nov 21 14:10:26 CET 2019
# Made by L-F-S
# At the University Of Trento, Italy

# un po' incasinato qui, devo decidere se processo solo l UI (tipoi in def
# feature, oppre se invece processo anche i data, come in tutte le altre tipo,
# e allora diventa un attimino più stile che fa tutto questa

#import os
import sys
import numpy as np
import re
import pandas as pd

#sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
#import filename_discrepancies


dataset="SmitsSA_2017" #sys.argv[1] #TODO switch test dataset
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" #TODO add nstepdir
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"


def feature(UI):
    if UI:
        return UI
    return "Cas9"

def subset_by_unknownSGB(data,UI):
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

def subset_by_SGB(data, UI):
    """ subset only data belonging to that SGB"""
    SGB=int(UI)
    return data[data["SGB ID"]==SGB]

def subset_by_species(data, UI):
    """subset data relative to that species
    accepted input (case insensitive):
            <Genus species>
            <any_taxonomical_grade>"""
    if UI:
        gen_spec=UI.split(" ")
        if len(gen_spec) == 2:
            genus=gen_spec[0][0].upper()+gen_spec[0][1:]
            species=gen_spec[1]
            line=genus+"_"+species
        elif len(gen_spec) == 1:
            line=gen_spec[0]
        else:
           raise Exception("ERROR: -species: \naccepted input (case INSENSITIVE):\n\t<Genus species>\n\t<any_taxonomical_grade>")

        return data[data["Estimated taxonomy"].str.contains(line, flags=re.IGNORECASE)]
    return data

def subste_by_active(data, UI):
    if UI:
        if UI.lower()=="n":
            return data
    return data.dropna(how="any")

def subset_by_length(data, UI):
    if UI:
        nmin=int(UI.split(" ")[0])
        nmax=int(UI.split(" ")[1])
        sorted_feature_counts=data["Seq"].str.count("").sort_values()
        length_intreval=sorted_feature_counts[sorted_feature_counts>=nmin]
        length_intreval=length_intreval[length_intreval<=nmax]
        subset=data.loc[length_intreval.index]
        return subset
    return data

def subste_by_genome(data, UI):
    return data[data.genome==UI]

def subset_by_type():
    return "TODO"

def plot_distributions(data):
    return "TODO"
