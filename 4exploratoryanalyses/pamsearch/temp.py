#!/bin/bash
# Mon Jul 29 10:24:14 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
sys.path.insert(0, '/home/lorenzo.signorini/utils/')
import filename_discrepancies


dataset="ZellerG_2014" #sys.argv[1] #TODO switch test dataset
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" #TODO add nstepdir
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"

oridir="/scratchCM/tmp_projects/epasolli_darkmatter/allcontigs/SmitsSA_2017/metabat/genomes_comp50_cont05/prokka/SmitsSA_2017__TZ_10742_megahit__bin.3"
filename=oridir+"/SmitsSA_2017__TZ_10742_megahit__bin.3.gff"
cmd="grep '^>' "+filename
contigs=subprocess.check_output(cmd, shell=True, universal_newlines=True).strip("\n")
print(type(contigs))
for unpolished_contig in contigs.split("\n"):
    print(unpolished_contig.lstrip(">"))
