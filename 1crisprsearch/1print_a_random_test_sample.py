# -*- coding: utf-8 -*-
"""
10/06/2019

@author L-F-S
at the University of Trento, Italy

Prints one random sample from one random study.
"""

import random
import os

gnms="/scratchCM/tmp_projects/epasolli_darkmatter/allcontigs/ALLreconstructedgenomes"
dirs=os.listdir(gnms)
os.chdir(gnms)
#%%

rand_fold=(dirs[random.randint(0,len(dirs)-1)])
print(rand_fold)
#print(os.listdir())
os.chdir(gnms+"/"+rand_fold)
bins=os.listdir()
rand_bin=(bins[random.randint(0,len(bins)-1)])
print(rand_bin)
rand_study=rand_bin.split("__")[0]
rand_sample=rand_bin.split("__")[1]
print("study: ",rand_study,"\nsample: ", rand_sample)
