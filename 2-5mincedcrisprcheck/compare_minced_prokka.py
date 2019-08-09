#!/bin/bash
# 21/06/2019
# Made by L-F-S
# At the University Of Trento, Italy

#compares the output of my minced vs prokkas' minced

#if they are equal, outputs nothing,
#else: outputs the difference as follows

# filename
#  my_minced: +++++
# <  roba>
#  prokka: ++++++
# < altraroba>

import sys
import os
import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt

outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2-5mincedcrisprcheck"
my_minced_dir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/1crisprsearch/out" #TODO leva test metti out
prokka_dir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno" #TODO leva test metti crisprcasanno
dataset = sys.argv[1] #TODO leva dataset zellerg (test)
mkdir="mkdir "+outdir+"/"+dataset
os.system(mkdir)

def compare_line_by_line(filename, my_processed_minced, prokka_processed_minced):
    print("comparing minced and prokka")
    
    os.system("diff "+my_processed_minced+" "+prokka_processed_minced+">"+outdir+"/"+dataset+"/"+filename+".dff") 
    return
def retrieve_prokka(filename):
    print("retrieving prokka annotation")
    cmd="less "+filename+" | grep minced >../justminced/"+dataset+"/"+filename+".minced"
    os.system(cmd)
    return prokka_dir.rstrip("/crisprcasanno")+"/justminced/"+dataset+"/"+filename+".minced"#TODO sostituisci/leva testdir

os.chdir(prokka_dir)
for filename in os.listdir():
    if filename.startswith(dataset):
        print("processing ",filename)
        processed_prokka=retrieve_prokka(filename)
        my_minced_file=my_minced_dir+"/"+filename.rstrip(".crisprcas.gff")+".fa.minced.out.gff"
        compare_line_by_line(filename.rstrip(".crisprcas.gff"),my_minced_file,processed_prokka)
