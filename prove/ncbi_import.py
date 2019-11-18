# Thu Nov 7 11:31:54 CET 2019
# Made by L-F-S
# At the University Of Trento, Italy

import os
import sys
import numpy as np
import pandas as pd
from Bio import Entrez


sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies


dataset="SmitsSA_2017" #sys.argv[1] #TODO switch test dataset
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" #TODO add nstepdir
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"

Entrez.email="my_mail"
handle=Entrez.einfo()
res=Entrez.read(handle)
print("as list:")
print(res['DbList'])

res=Entrez.read(Entrez.einfo(db= "sra"))

print(res["DbInfo"]["Description"])
