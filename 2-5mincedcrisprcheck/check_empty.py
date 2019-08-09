#!/bin/bash
# 06/06/2019
# Made by L-F-S
# At the University Of Trento, Italy

import os
import numpy as np
import pandas as pd
import sys
testdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2-5mincedcrisprcheck/ZellerG_2014"
dataset_dir=sys.argv[1]
directory="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2-5mincedcrisprcheck/"+dataset_dir
os.chdir(directory)
for filex in os.listdir():
    f=open(filex)
    file1nodes=[]
    effettivamentediversi=False
    for line in f.readlines():
        if line.startswith("<"):  #denotes line in file 1
            file1nodes.append(line.split()[1])
        if line.startswith(">"): #denotes line  in file 2
            file2node=line.split()[1]
            if not file2node in file1nodes:
                effettivamentediversi=True
    if not effettivamentediversi:
        cmd= "rm "+filex
        os.system(cmd)
        print("not really different, removing file ", filex)
