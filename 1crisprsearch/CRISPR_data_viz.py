#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 18:26:11 2019

@author: L-F-S
@ University of Trento, Italy
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Extract summary features from CRISPR data
load_data=False


if load_data==False:
    data_path="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/1crisprsearch/out"
    os.chdir(data_path)
    def add_ncrispr_from_pilercr(bins, filename):
        f=open(filename)
        for line in f.readlines():
            if "putative CRISPR" in line:
                bins[filename.rstrip(".fa.pilercr.out")][1]= int(line.split(" ")[1])       
        f.close()
        return bins

    def add_ncrispr_from_minced(bins, filename):
        f=open(filename)
        n=0
        for line in f.readlines():
            if line.startswith("NODE"):
                n+=1
        bins[filename.rstrip(".fa.minced.out.gff")][2]= n
        f.close()
        return bins

    # times are calculated differently so they're not relevant, ma vabbe
    def add_pilercr_time(bins, filename):
        f=open(filename)
        bins[filename.rstrip(".fa.pilercr.time")][3]= float(f.readline().split(" ")[1].rstrip("\n"))        
        f.close()
        return bins

    def add_minced_time(bins,filename):
        f=open(filename)
        n=0
        for line in f.readlines():
            if line.startswith("Sequence"):
                n+=1
            if line.startswith("Time"):
                bins[filename.rstrip(".fa.minced.out")][4] = float(line.split(" ")[4])/1000
        f.close()
        return bins
        
    #  CREA TABELLONA GIABIGA

    bins={}

    for filename in os.listdir():
        print(filename)
        if "pilercr" in filename:
            if not "time" in filename: # it's the pilercr output file 
                if not filename.rstrip(".fa.pilercr.out") in bins.keys(): # add binname to dictionary  
                    bins[filename.rstrip(".fa.pilercr.out")]=[filename.rstrip(".fa.pilercr.out"),0,0,0,0]
                bins=add_ncrispr_from_pilercr(bins, filename)
                    
            else:                  # it's the pilercr execution time
                if not filename.rstrip(".fa.pilercr.time") in bins.keys(): # add binname to dictionary  
                    bins[filename.rstrip(".fa.pilercr.time")]=[filename.rstrip(".fa.pilercr.time"),0,0,0,0]
    #            bins=add_pilercr_time(bins, filename)
        
        else:   
            if ".gff" in filename: # it's the minced .gff file
                if not filename.rstrip(".fa.minced.out.gff") in bins.keys():  # add binname to dictionary #TODO : leva il .txt dalle estensioni file finali
                    bins[filename.rstrip(".fa.minced.out.gff")]=[filename.rstrip(".fa.minced.out.gff"),0,0,0,0]
                bins=add_ncrispr_from_minced(bins, filename)
            
            else:                  # it's the minced output file
                if not filename.rstrip(".fa.minced.out") in bins.keys():  # add binname to dictionary #TODO : leva il .txt dalle estensioni file finali
                    bins[filename.rstrip(".fa.minced.out")]=[filename.rstrip(".fa.minced.out"),0,0,0,0]
     #           bins=add_minced_time(bins, filename)

    # convertila a pd.DataFrame
    df = pd.DataFrame.from_dict(bins, orient='index', columns = ["name","nCRISPR_pilercr","nCRISPR_minced","t_pilercr","t_minced"])
    df.to_csv("../CRISPRhisto.csv")


else: #load_data==False
    print("loading data from file")
    os.chdir("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/1crisprsearch/")
    df=pd.read_csv("CRISPRhisto.csv")

# total crispr found
l0= str(df.shape[0]) + "bins (=reconstructed genomes) analyzed"
l1 = str(sum(df["nCRISPR_minced"])) + " CRISPR sequences found with minced"
l2 = str(sum(df["nCRISPR_pilercr"])) + " CRISPR sequences founf with pilercr"

f=open("../totalcrispr", "w")
f.write(l0+"\n"+l1+"\n"+l2+"\n")
f.close()

# PLOT

plt.figure()
plt.hist(df["nCRISPR_minced"], color="orange", alpha=0.4, label="minced")
plt.hist(list(df["nCRISPR_pilercr"]), color="teal", alpha=0.4, label="pilercr")
plt.legend()
plt.savefig("../CRISPRhist.png")
