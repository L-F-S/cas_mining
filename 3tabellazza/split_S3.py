#!/bin/bash
# 06/06/2019
# Made by L-F-S
# At the University Of Trento, Italy
# Routine to split $data/S3Segata.py in
# N different .csv , where N is the nubmer of studies.
# this will be useful later to parallelize table creation

import os
import sys
import numpy as np
import pandas as pd


dataset="ZellerG_2014" #sys.argv[1]
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
CRISPRdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/1crisprsearch/out"
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"

bigtable=pd.read_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/S3Segata.csv")

for study, data_study in bigtable.groupby(bigtable.Study):
    print(study,data_study.shape)
    data_study.to_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/S3Segata_split/S3_subs_"+study)
