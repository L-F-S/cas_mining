# Thu May 7 15:06:09 CEST 2020
# Made by L-F-S
# At the University Of Trento, Italy

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/9casminer_v2/IO/')
sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/9casminer_v2/')
import filename_discrepancies
import locus


dataset="SmitsSA_2017" #sys.argv[1] #TODO switch test dataset
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" #TODO add nstepdir
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"

""" takes as in
